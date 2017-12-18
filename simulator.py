from IndividualSimulator.code import config
from IndividualSimulator.code import hosts
from IndividualSimulator.code import outputdata
from IndividualSimulator.code.eventhandling import EventHandler
from IndividualSimulator.code.interventionhandling import InterventionHandler
from IndividualSimulator.code.ratehandling import RateHandler
from IndividualSimulator.code.ratestructures.ratetree import RateTree
import argparse
import copy
import inspect
import numpy as np
import time as time_mod
import raster_tools


__version__ = "0.0.10"


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

        # Setup function to get next state from current state
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

        # Read in hosts
        self.params['init_hosts'], self.params['init_cells'] = hosts.read_host_files(
                self.params['HostPosFile'].split(","), self.params['InitCondFile'].split(","),
                self.params['RegionFile'], states, sim_type=self.params['SimulationType'])

        self.params['nhosts'] = len(self.params['init_hosts'])
        if self.params['init_cells'] is not None:
            self.params['ncells'] = len(self.params['init_cells'])

        # Kernel setup
        if self.params['KernelType'] == "EXPONENTIAL":
            self.params['kernel'] = kernel_exp(self.params['KernelScale'])
        elif self.params['KernelType'] == "NONSPATIAL":
            self.params['kernel'] = kernel_nonspatial()
        elif self.params['KernelType'] == "RASTER":
            self.params['kernel'] = raster_tools.RasterData.from_file(
                self.params['KernelFile']).array
        else:
            raise ValueError("Unrecognised KernelType!")

        # Setup Virtual Sporulation
        if self.params['SimulationType'] == "RASTER":
            if self.params['VirtualSporulationStart'] is None:

                array_size = self.params['kernel'].shape
                self.params['coupled_positions'] = [
                    (x, y) for y in range(-array_size[1], array_size[1]+1)
                    for x in range(-array_size[0], array_size[0]+1)]
                
                self.params['coupled_kernel'] = self.params['kernel']

            else:
                start = self.params['VirtualSporulationStart']
                self.params['coupled_positions'] = [
                    (x, y) for y in range(1-start, start)
                    for x in range(1-start, start)]

                self.params['coupled_kernel'] = self.params['kernel'][0:start, 0:start]


                vs_kernel = np.copy(self.params['kernel'])
                vs_kernel[0:start, 0:start] = 0
                spore_prob = np.sum(vs_kernel)
                vs_kernel = vs_kernel.flatten() / spore_prob
                self.params['vs_kernel'] = RateTree(len(vs_kernel))
                for i, kernel_val in enumerate(vs_kernel):
                    self.params['vs_kernel'].insert_rate(i, kernel_val)


                self.params['spore_rate'] = self.params['InfRate'] * spore_prob

        # self.params['init_region_summary'] = [{key: 0 for key in states + ["Culled"]}
        #                                       for _ in range(self.params['NRegions'])]

        self.rate_handler = RateHandler(self)

        # Setup initial rates
        if self.params['CacheKernel'] is True and self.params['SimulationType'] == "INDIVIDUAL":
            distances, kernel_vals = calc_dist_kernel(self.params['init_hosts'],
                                                      self.params['kernel'])
            self.params['kernel_vals'] = kernel_vals
            self.params['distances'] = distances

        self.event_handler = EventHandler(self, self.rate_handler)

        self.params['region_map'] = {key: [] for key in range(self.params['NRegions'])}
        if self.params['init_cells'] is None:
            self.params['init_inf_rates'] = np.zeros(self.params['nhosts'])
        else:
            self.params['init_inf_rates'] = np.zeros(self.params['ncells'])
            if self.params['VirtualSporulationStart'] is not None:
                self.params['init_spore_rates'] = np.zeros(self.params['ncells'])
        self.params['init_adv_rates'] = np.zeros(self.params['nhosts'])

        if self.params['SimulationType'] == "INDIVIDUAL":

            for i in range(self.params['nhosts']):
                current_state = self.params['init_hosts'][i].state
                region = self.params['init_hosts'][i].reg
                # self.params['init_region_summary'][region][current_state] += 1
                self.params['region_map'][region].append(i)
                if current_state in "ECDI":
                    self.params['init_adv_rates'][i] = self.params[current_state + 'AdvRate']
                    if current_state in "CI":
                        for j in range(self.params['nhosts']):
                            if self.params['init_hosts'][j].state == "S":
                                self.params['init_inf_rates'][j] += self.event_handler.kernel(j, i)

        elif self.params['SimulationType'] == "RASTER":

            self.params['cell_map'] = {}

            for cell in self.params['init_cells']:
                self.params['cell_map'][cell.cell_position] = cell.cell_id

            for cell in self.params['init_cells']:
                for host in cell.hosts:
                    current_state = host.state
                    region = host.reg
                    self.params['region_map'][region].append(host.host_id)
                    if current_state in "ECDI":
                        self.params['init_adv_rates'][host.host_id] = self.params[
                            current_state + 'AdvRate']
                if (cell.states["C"] + cell.states["I"]) > 0:
                    for cell2_rel_pos in self.params['coupled_positions']:
                        cell2_pos = tuple(item1 + item2 for item1, item2
                                          in zip(cell.cell_position, cell2_rel_pos))
                        cell2_id = self.params['cell_map'].get(cell2_pos, None)
                        if cell2_id is None:
                            continue
                        cell2 = self.params['init_cells'][cell2_id]
                        self.params['init_inf_rates'][cell2_id] += cell2.states["S"] * (
                            (cell.states["C"] + cell.states["I"]) *
                            self.event_handler.kernel(cell2_rel_pos) / 100)

                    if self.params['VirtualSporulationStart'] is not None:
                        self.params['init_spore_rates'][cell.cell_id] = (cell.states["C"] +
                                                                         cell.states["I"])

        else:
            raise ValueError("Unrecognised SimulationType!")

        self.rate_factor = [self.params['InfRate'], 1]

        if self.params['VirtualSporulationStart'] is not None:
            self.rate_factor.append(self.params['spore_rate'])

        # Intervention setup
        self.intervention_handler = InterventionHandler(self)

        end_time = time_mod.time()

        print(len(self.params['init_hosts']),len(self.params['init_cells']))

        if silent is False:
            print("Initial setup complete.  "
                  "Time taken: {0:.3f} seconds.".format(end_time - start_time))

    def run_epidemic(self, iteration=0, silent=False):
        start_time = time_mod.time()

        # Setup run parameters to keep track of iteration
        self.run_params = {}
        self.run_params['all_events'] = []
        # self.run_params['region_summary'] = copy.deepcopy(self.params['init_region_summary'])
        # self.run_params['summary_dump'] = []

        # self.all_hosts = copy.deepcopy(self.params['init_hosts'])
        # self.all_cells = copy.deepcopy(self.params['init_cells'])
        self.all_hosts = self.params['init_hosts']
        self.all_cells = self.params['init_cells']

        print("Copied hosts and cells", flush=True)

        # Zero all rates
        self.rate_handler.zero_rates()

        # Initialise rates from setup
        if self.params['SimulationType'] == "INDIVIDUAL":
            for i in range(self.params['nhosts']):
                self.rate_handler.insert_rate(i, self.params['init_inf_rates'][i], "Infection")
                self.rate_handler.insert_rate(i, self.params['init_adv_rates'][i], "Advance")

        elif self.params['SimulationType'] == "RASTER":
            for i in range(self.params['nhosts']):
                self.rate_handler.insert_rate(i, self.params['init_adv_rates'][i], "Advance")

            for j in range(self.params['ncells']):
                self.rate_handler.insert_rate(j, self.params['init_inf_rates'][j], "Infection")
                if self.params['VirtualSporulationStart'] is not None:
                    self.rate_handler.insert_rate(
                        j, self.params['init_spore_rates'][j], "Sporulation")

        else:
            raise ValueError("Unrecognised SimulationType!")

        self.time = 0
        # self.run_params['summary_dump'].append(
        #     (self.time, copy.deepcopy(self.params['init_region_summary'])))
        # if self.params['SummaryOutputFreq'] != 0:
        #     nextSummaryDumpTime = self.time + self.params['SummaryOutputFreq']
        # else:
        #     nextSummaryDumpTime = np.inf

        # Initialise interventions
        self.intervention_handler.initialise_rates(self.all_hosts)

        # Set time until first intervention
        nextInterventionTime = self.intervention_handler.next_intervention_time

        # Run gillespie loop
        while True:
            # Find next event from event handler
            totRate, event_type, hostID = self.rate_handler.get_next_event()
            if event_type is None:
                nextTime = np.inf
            else:
                nextTime = self.time + (-1.0/totRate)*np.log(np.random.random_sample())

            if nextTime >= nextInterventionTime:
                if nextInterventionTime > self.params['FinalTime']:
                    break
                self.time = nextInterventionTime
                # carry out intervention update
                self.intervention_handler.update(self.all_hosts, self.time)
                nextInterventionTime = self.intervention_handler.next_intervention_time
                # else:
                #     if nextSummaryDumpTime >= self.params['FinalTime']:
                #         break
                #     # Dump summary data
                #     self.time = nextSummaryDumpTime
                #     self.run_params['summary_dump'].append(
                #         (nextSummaryDumpTime, copy.deepcopy(self.run_params['region_summary'])))
                #     nextSummaryDumpTime += self.params['SummaryOutputFreq']
            else:
                if nextTime > self.params['FinalTime']:
                    break
                # Carry out event
                self.time = nextTime
                self.event_handler.do_event(event_type, hostID, self.all_hosts, self.all_cells)
                if self.params['UpdateOnAllEvents'] is True:
                    self.intervention_handler.update_on_event(self.all_hosts, self.time)
                    nextInterventionTime = self.intervention_handler.next_intervention_time

        # self.run_params['summary_dump'].append(
        #     (nextSummaryDumpTime, copy.deepcopy(self.run_params['region_summary'])))

        end_time = time_mod.time()

        if silent is False:
            print("Run {0} of {1} complete.  ".format(iteration+1, self.params['NIterations']) +
                  "Time taken: {0:.3f} seconds.".format(end_time - start_time), end="\n")

        print("Total number of cells infected: {0}".format(np.sum([1 for x in self.all_cells if x.states["I"] > 0])))

        return (self.all_hosts, self.all_cells, self.run_params)

    def output_run_data(self, all_hosts, all_cells, run_params, iteration=0):
        run_data = outputdata.output_all_run_data(self, all_hosts, all_cells, run_params, iteration)
        return run_data


def run_epidemics(params):
    all_data = []
    if params['SaveSetup']:
        run_sim = Simulator(params)
        run_sim.setup()
    for iteration in range(params['NIterations']):
        if not params['SaveSetup']:
            run_sim = Simulator(params)
            run_sim.setup()
        all_hosts, all_cells, run_params = run_sim.run_epidemic(iteration)
        all_data.append(
            run_sim.output_run_data(all_hosts, all_cells, run_params, iteration=iteration))

    return all_data


def main(configFile="config.ini", keyFile=False, defaultConfig=None, params_options=None):
    frame = inspect.stack()[1]
    modu = inspect.getmodule(frame[0])

    if keyFile is True:
        config.write_keyfile()
        print("KeyFile generated.")

    if defaultConfig is not None:
        config.write_default_config(defaultConfig)
        print("Default config file generated.")

    if keyFile is False and defaultConfig is None:
        params = config.read_config_file(filename=configFile)
        if params_options is not None:
            for key, value in params_options.items():
                # TODO Should check that these are valid additions
                params[key] = value

        config.check_params_valid(params)

        params['call_params'] = copy.deepcopy(params)
        params['call_config_file'] = configFile
        params['call_module'] = str(modu)
        params['call_time'] = time_mod.strftime("%a, %d %b %Y %H:%M:%S", time_mod.localtime())
        params['call_version'] = __version__
        all_data = run_epidemics(params)

        outputdata.output_log_file(params)

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
