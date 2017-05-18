import config
import argparse
import eventhandling
import outputdata
import copy
import time as time_mod
import numpy as np
from host import Host
from ratesum import RateSum
from rateinterval import RateInterval


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

        self.params['init_region_summary'] = [{key: 0 for key in states}
                                              for _ in range(self.params['NRegions'])]

        if self.params['RateStructure-Infection'] == "ratesum":
            self.inf_rates = RateSum(self.params['nhosts'])
        elif self.params['RateStructure-Infection'] == "rateinterval":
            self.inf_rates = RateInterval(self.params['nhosts'])
        else:
            raise ValueError("Invalid rate structure - infection events!")

        if self.params['RateStructure-Advance'] == "ratesum":
            self.adv_rates = RateSum(self.params['nhosts'])
        elif self.params['RateStructure-Advance'] == "rateinterval":
            self.adv_rates = RateInterval(self.params['nhosts'])
        else:
            raise ValueError("Invalid rate structure - advance events!")

        self.all_rates = [self.inf_rates, self.adv_rates]

        # Setup initial rates
        if self.params['CacheKernel'] is True:
            distances, kernel_vals = calc_dist_kernel(self.params['init_hosts'],
                                                      self.params['kernel'])
            self.params['kernel_vals'] = kernel_vals
            self.params['distances'] = distances
            self.do_event = eventhandling.do_event_cached

        else:
            self.do_event = eventhandling.do_event_uncached

        self.params['init_inf_rates'] = np.zeros(self.params['nhosts'])
        self.params['init_adv_rates'] = np.zeros(self.params['nhosts'])
        for i in range(self.params['nhosts']):
            current_state = self.params['init_hosts'][i].state
            region = self.params['init_hosts'][i].reg
            self.params['init_region_summary'][region][current_state] += 1
            if current_state in "ECDI":
                self.params['init_adv_rates'][i] = self.params[current_state + 'AdvRate']
                if current_state in "CI":
                    for j in range(self.params['nhosts']):
                        if self.params['init_hosts'][j].state == "S":
                            if self.params['CacheKernel'] is True:
                                self.params['init_inf_rates'][j] += self.params[
                                    'kernel_vals'][j, i]
                            else:
                                self.params['init_inf_rates'][j] += self.params['kernel'](
                                    np.linalg.norm([
                                        self.params['init_hosts'][i].x -
                                        self.params['init_hosts'][j].x,
                                        self.params['init_hosts'][i].y -
                                        self.params['init_hosts'][j].y]))

        end_time = time_mod.time()

        if silent is False:
            print("Initial setup complete.  "
                  "Time taken: {0:.3f} seconds.".format(end_time - start_time))

    def run_epidemic(self, iteration=0):
        start_time = time_mod.time()

        run_params = {}
        run_params['all_events'] = []
        run_params['region_summary'] = copy.deepcopy(self.params['init_region_summary'])
        run_params['summary_dump'] = []

        all_hosts = copy.deepcopy(self.params['init_hosts'])

        for rates in self.all_rates:
            rates.zero_rates()

        for i in range(self.params['nhosts']):
            self.inf_rates.insert_rate(i, self.params['init_inf_rates'][i])
            self.adv_rates.insert_rate(i, self.params['init_adv_rates'][i])

        time2 = time_mod.time()

        time = 0
        run_params['summary_dump'].append(
            (time, copy.deepcopy(self.params['init_region_summary'])))
        if self.params['SummaryOutputFreq'] != 0:
            nextSummaryDumpTime = time + self.params['SummaryOutputFreq']
        else:
            nextSummaryDumpTime = self.params['FinalTime']

        # Run gillespie loop
        while True:
            self.inf_rates.full_resum()  # As errors accumulate in total rate
            self.adv_rates.full_resum()  # As errors accumulate in total rate
            totRates = [self.params['InfRate']*self.inf_rates.get_total_rate(),
                        self.adv_rates.get_total_rate()]
            totRate = np.sum(totRates)
            if totRate <= 0:
                break

            nextTime = (-1.0/totRate)*np.log(np.random.random_sample())
            select_rate = np.random.random_sample()*totRate

            time += nextTime

            while time >= nextSummaryDumpTime:
                run_params['summary_dump'].append((nextSummaryDumpTime,
                                                  copy.deepcopy(run_params['region_summary'])))
                nextSummaryDumpTime += self.params['SummaryOutputFreq']

            if time > self.params['FinalTime']:
                break

            if select_rate < totRates[0]:
                eventID = self.inf_rates.select_event(select_rate/self.params['InfRate'])
                self.do_event(eventID, all_hosts, self.all_rates, self.params, run_params, time)
            elif select_rate < np.sum(totRates):
                eventID = self.adv_rates.select_event(
                    select_rate - self.params['InfRate']*self.inf_rates.get_total_rate())
                self.do_event(eventID, all_hosts, self.all_rates, self.params, run_params, time)

        run_params['summary_dump'].append((nextSummaryDumpTime,
                                          copy.deepcopy(run_params['region_summary'])))

        end_time = time_mod.time()

        print("Run {0} of {1} complete.  ".format(iteration+1, self.params['NIterations']) +
              "Time taken: {0:.3f} seconds.".format(end_time - start_time))

        return (all_hosts, run_params)

    def output_run_data(self, all_hosts, run_params, iteration=0):
        outputdata.output_data_hosts(all_hosts, self.params, iteration=iteration)
        outputdata.output_data_events(run_params, iteration=iteration)
        if self.params['SummaryOutputFreq'] != 0:
            outputdata.output_data_summary(self.params, run_params, iteration=iteration)


def run_epidemics(params):
    run_sim = Simulator(params)
    run_sim.setup()
    for iteration in range(params['NIterations']):
        all_hosts, run_params = run_sim.run_epidemic(iteration)
        run_sim.output_run_data(all_hosts, run_params, iteration=iteration)


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

    if args.keyFile is True:
        config.write_keyfile()
        print("KeyFile generated.")

    if args.defaultConfig is not None:
        config.write_default_config(args.defaultConfig)
        print("Default config file generated.")

    if args.keyFile is False and args.defaultConfig is None:
        params = config.read_config_file(filename=args.configFile)
        run_epidemics(params)
