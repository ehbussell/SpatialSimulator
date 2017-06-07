import config
import argparse
import eventhandling
from eventhandling import EventHandler
import outputdata
import copy
import importlib
import time as time_mod
import numpy as np
from host import Host
from ratehandling import RateHandler


def kernel_exp(kernel_param):

    def kernel(dist):
        if dist > 0:
            return np.exp(-kernel_param*dist)
        else:
            return 0

    return kernel


def kernel_nonspatial():

    def kernel(dist):
        return 1.0

    return kernel


def calc_dist_kernel(hosts, kernel):
    nhosts = len(hosts)
    distances = np.zeros((nhosts, nhosts))
    kernel_vals = np.zeros((nhosts, nhosts))

    for i in range(nhosts):
        for j in range(i):
            dist = np.linalg.norm([hosts[i].x-hosts[j].x, hosts[i].y-hosts[j].y])
            distances[i, j] = dist
            distances[j, i] = dist

            kval = kernel(dist)
            kernel_vals[i, j] = kval
            kernel_vals[j, i] = kval

    return (distances, kernel_vals)


class Simulator:

    def __init__(self, params=None, config_file=None):
        if params is not None:
            self.params = copy.deepcopy(params)
        elif config_file is not None:
            params = config.read_config_file(filename=config_file)
            self.params = copy.deepcopy(params)
        else:
            raise ValueError("Must specify either parameter dictionary or configuration file!")

    def setup(self, silent=False):
        """Run all static setup before epidemic simulations."""

        start_time = time_mod.time()

        # Read in hosts
        self.params['init_hosts'] = config.read_hosts(self.params['HostFile'])
        self.params['nhosts'] = len(self.params['init_hosts'])

        # Setup
        states = list(self.params['Model'])

        def next_state_func(current_state):
            states_iter = iter(states)
            state = next(states_iter, None)
            while current_state != state:
                state = next(states_iter, None)
                if state is None:
                    break
            return next(states_iter, None)

        self.params['next_state'] = next_state_func

        if self.params['KernelType'] == "EXPONENTIAL":
            self.params['kernel'] = kernel_exp(self.params['KernelScale'])
        elif self.params['KernelType'] == "NONSPATIAL":
            self.params['kernel'] = kernel_nonspatial()
        else:
            raise ValueError("Unrecognised KernelType!")

        self.params['init_region_summary'] = [{key: 0 for key in states + ["Culled"]}
                                              for _ in range(self.params['NRegions'])]

        self.rate_handler = RateHandler(self)

        # Setup initial rates
        if self.params['CacheKernel'] is True:
            distances, kernel_vals = calc_dist_kernel(self.params['init_hosts'],
                                                      self.params['kernel'])
            self.params['kernel_vals'] = kernel_vals
            self.params['distances'] = distances

        self.event_handler = EventHandler(self, self.rate_handler)

        self.params['region_map'] = {key: [] for key in range(self.params['NRegions'])}
        self.params['init_inf_rates'] = np.zeros(self.params['nhosts'])
        self.params['init_adv_rates'] = np.zeros(self.params['nhosts'])
        for i in range(self.params['nhosts']):
            current_state = self.params['init_hosts'][i].state
            region = self.params['init_hosts'][i].reg
            self.params['init_region_summary'][region][current_state] += 1
            self.params['region_map'][region].append(i)
            if current_state in "ECDI":
                self.params['init_adv_rates'][i] = self.params[current_state + 'AdvRate']
                if current_state in "CI":
                    for j in range(self.params['nhosts']):
                        if self.params['init_hosts'][j].state == "S":
                            self.params['init_inf_rates'][j] += self.event_handler.kernel(j, i)

        end_time = time_mod.time()

        if silent is False:
            print("Initial setup complete.  "
                  "Time taken: {0:.3f} seconds.".format(end_time - start_time))

    def run_epidemic(self, iteration=0, silent=False):
        start_time = time_mod.time()

        self.run_params = {}
        self.run_params['all_events'] = []
        self.run_params['region_summary'] = copy.deepcopy(self.params['init_region_summary'])
        self.run_params['summary_dump'] = []

        self.all_hosts = copy.deepcopy(self.params['init_hosts'])

        self.rate_factor = [self.params['InfRate'], 1]
        self.rate_handler.zero_rates()

        for i in range(self.params['nhosts']):
            self.rate_handler.insert_rate(i, self.params['init_inf_rates'][i], "Infection")
            self.rate_handler.insert_rate(i, self.params['init_adv_rates'][i], "Advance")

        self.time = 0
        self.run_params['summary_dump'].append(
            (self.time, copy.deepcopy(self.params['init_region_summary'])))
        if self.params['SummaryOutputFreq'] != 0:
            nextSummaryDumpTime = self.time + self.params['SummaryOutputFreq']
        else:
            nextSummaryDumpTime = np.inf

        # Setup Interventions
        if self.params['InterventionScript'] is not None:
            # import and create intervention class
            intervention_module = importlib.import_module(self.params['InterventionScript'])
            self.params['intervention'] = intervention_module.create_interventions(self)
        else:
            self.params['intervention'] = None

        if self.params['intervention'] is None:
            nextInterventionTime = np.inf
        else:
            nextInterventionTime = self.time + self.params['intervention'].update_freq

        # Run gillespie loop
        while True:
            # Find next event from event handler
            totRate, event_type, hostID = self.rate_handler.get_next_event()
            if event_type is None:
                nextTime = np.inf
            else:
                nextTime = self.time + (-1.0/totRate)*np.log(np.random.random_sample())

            if nextTime >= nextSummaryDumpTime or nextTime >= nextInterventionTime:
                if nextSummaryDumpTime >= nextInterventionTime:
                    if nextInterventionTime > self.params['FinalTime']:
                        break
                    self.time = nextInterventionTime
                    # carry out intervention update
                    self.params['intervention'].update()
                    nextInterventionTime += self.params['intervention'].update_freq
                else:
                    if nextSummaryDumpTime >= self.params['FinalTime']:
                        break
                    # Dump summary data
                    self.time = nextSummaryDumpTime
                    self.run_params['summary_dump'].append(
                        (nextSummaryDumpTime, copy.deepcopy(self.run_params['region_summary'])))
                    nextSummaryDumpTime += self.params['SummaryOutputFreq']
            else:
                if nextTime > self.params['FinalTime']:
                    break
                # Carry out event
                self.time = nextTime
                self.event_handler.do_event(event_type, hostID, self.all_hosts)
                if self.params['UpdateOnAllEvents'] is True:
                    self.params['intervention'].update()

        self.run_params['summary_dump'].append(
            (nextSummaryDumpTime, copy.deepcopy(self.run_params['region_summary'])))

        end_time = time_mod.time()

        if silent is False:
            print("Run {0} of {1} complete.  ".format(iteration+1, self.params['NIterations']) +
                  "Time taken: {0:.3f} seconds.".format(end_time - start_time), end="\n")

        return (self.all_hosts, self.run_params)

    def output_run_data(self, all_hosts, run_params, iteration=0):
        run_data = outputdata.output_all_run_data(self, all_hosts, run_params, iteration)
        return run_data


def run_epidemics(params):
    all_data = []
    run_sim = Simulator(params)
    run_sim.setup()
    for iteration in range(params['NIterations']):
        all_hosts, run_params = run_sim.run_epidemic(iteration)
        all_data.append(run_sim.output_run_data(all_hosts, run_params, iteration=iteration))

    return all_data


def main(configFile="config.ini", keyFile=False, defaultConfig=None):
    if keyFile is True:
        config.write_keyfile()
        print("KeyFile generated.")

    if defaultConfig is not None:
        config.write_default_config(defaultConfig)
        print("Default config file generated.")

    if keyFile is False and defaultConfig is None:
        params = config.read_config_file(filename=configFile)
        all_data = run_epidemics(params)

        return all_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-c", "--configFile", default="config.ini",
                        help="Name of config file to use.")
    parser.add_argument("-k", "--keyFile", action="store_true",
                        help="Flag to generate keyfile.  "
                        "If present keyfile is created and program exits.")
    parser.add_argument("-d", "--defaultConfig", const="config.ini",
                        nargs="?", help="Flag to generate default config file."
                        "  If present file is created and program exits.  "
                        "Optionally specify filename.")
    args = parser.parse_args()

    main(**vars(args))
