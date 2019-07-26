import pdb
from IndividualSimulator.code import config
from IndividualSimulator.code import hosts
from IndividualSimulator.code import outputdata
from IndividualSimulator.code.eventhandling import EventHandler
from IndividualSimulator.code.interventionhandling import InterventionHandler
from IndividualSimulator.code.ratehandling import RateHandler
from IndividualSimulator.code.ratestructures.ratetree import RateTree
from IndividualSimulator.code.ratestructures.ratesum import RateSum
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
            dist = np.linalg.norm([hosts[i].xpos-hosts[j].xpos, hosts[i].ypos-hosts[j].ypos])
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
        init_hosts, init_cells, header = hosts.read_host_files(
            self.params['HostPosFile'].split(","), self.params['InitCondFile'].split(","),
            self.params['RegionFile'], states, sim_type=self.params['SimulationType'])
        self.params['init_hosts'] = init_hosts
        self.params['init_cells'] = init_cells
        self.params['header'] = header

        # Read in susceptibility and infectiousness
        hosts.read_sus_inf_files(self.params['init_cells'], self.params['header'],
                                 self.params['SusceptibilityFile'],
                                 self.params['InfectiousnessFile'],
                                 sim_type=self.params['SimulationType'])

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
                    (x, y) for y in range(-int(array_size[1]/2), int(array_size[1]/2)+1)
                    for x in range(-int(array_size[0]/2), int(array_size[0]/2)+1)]

                self.params['coupled_kernel'] = self.params['kernel']

            else:
                start = self.params['VirtualSporulationStart']
                self.params['coupled_positions'] = [
                    (x, y) for y in range(1-start, start)
                    for x in range(1-start, start)]

                centre = [int(x/2) for x in self.params['kernel'].shape]

                self.params['coupled_kernel'] = self.params['kernel'][
                    (centre[0]+1-start):(centre[0]+start), (centre[1]+1-start):(centre[1]+start)]

                vs_kernel = np.copy(self.params['kernel'])
                vs_kernel[(centre[0]+1-start):(centre[0]+start),
                          (centre[1]+1-start):(centre[1]+start)] = 0
                spore_prob = np.sum(vs_kernel)
                vs_kernel = vs_kernel.flatten() / spore_prob
                self.params['vs_kernel'] = RateTree(len(vs_kernel))
                for i, kernel_val in enumerate(vs_kernel):
                    self.params['vs_kernel'].insert_rate(i, kernel_val)

                self.params['spore_rate'] = self.params['InfRate'] * spore_prob

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
                        self.params['init_inf_rates'][cell2_id] += (
                            cell2.susceptibility * cell2.states["S"] *
                            (cell.states["C"] + cell.states["I"]) * cell.infectiousness *
                            self.event_handler.kernel(cell2_rel_pos)) / self.params['MaxHosts']

                    if self.params['VirtualSporulationStart'] is not None:
                        self.params['init_spore_rates'][cell.cell_id] = (
                            cell.states["C"] + cell.states["I"]) * cell.infectiousness

        else:
            raise ValueError("Unrecognised SimulationType!")

        self.rate_factor = [self.params['InfRate'], 1]

        if self.params['VirtualSporulationStart'] is not None:
            self.rate_factor.append(self.params['spore_rate'])

        # Intervention setup
        self.intervention_handler = InterventionHandler(self)

        end_time = time_mod.time()

        if not silent:
            if self.params['init_hosts'] is not None:
                print("Num init hosts: ", len(self.params['init_hosts']))
            if self.params['init_cells'] is not None:
                print("Num init cells: ", len(self.params['init_cells']))
            print("Initial setup complete.  "
                  "Time taken: {0:.3f} seconds.".format(end_time - start_time))

    def initialise(self, silent=False):
        # Setup run parameters to keep track of iteration
        self.run_params = {}
        self.run_params['all_events'] = []
        # self.run_params['region_summary'] = copy.deepcopy(self.params['init_region_summary'])
        # self.run_params['summary_dump'] = []

        # self.all_hosts = copy.deepcopy(self.params['init_hosts'])
        # self.all_cells = copy.deepcopy(self.params['init_cells'])
        self.all_hosts = self.params['init_hosts']
        self.all_cells = self.params['init_cells']

        if not silent:
            print("Copied hosts and cells", flush=True)

        # Zero all rates
        self.rate_handler.zero_rates()

        # Initialise rates from setup in bulk
        if self.params['SimulationType'] == "INDIVIDUAL":
            self.rate_handler.bulk_insert(self.params['init_inf_rates'], "Infection")
            self.rate_handler.bulk_insert(self.params['init_adv_rates'], "Advance")

        elif self.params['SimulationType'] == "RASTER":
            self.rate_handler.bulk_insert(self.params['init_adv_rates'], "Advance")

            self.rate_handler.bulk_insert(self.params['init_inf_rates'], "Infection")
            if self.params['VirtualSporulationStart'] is not None:
                self.rate_handler.bulk_insert(self.params['init_spore_rates'], "Sporulation")


        else:
            raise ValueError("Unrecognised SimulationType!")

        self.time = 0

        # Initialise interventions
        self.intervention_handler.initialise_rates(self.all_hosts, self.all_cells)


    def run_epidemic(self, iteration=0, silent=False):
        start_time = time_mod.time()

        # Set time until first intervention
        next_intervention_time = self.intervention_handler.next_intervention_time
        if self.params['RasterOutputFreq'] != 0:
            outputdata.output_raster_data(self, time=self.time, iteration=iteration)
            nextRasterDumpTime = self.time + self.params['RasterOutputFreq']
        else:
            nextRasterDumpTime = np.inf

        # Run gillespie loop
        while True:
            # Find next event from event handler
            totRate, event_type, hostID = self.rate_handler.get_next_event()
            if event_type is None:
                nextTime = np.inf
            else:
                nextTime = self.time + (-1.0/totRate)*np.log(np.random.random_sample())

            while np.minimum(nextTime, next_intervention_time) >= nextRasterDumpTime:
                if nextRasterDumpTime > self.params['FinalTime']:
                    break
                self.time = nextRasterDumpTime
                outputdata.output_raster_data(self, time=self.time, iteration=iteration)
                nextRasterDumpTime += self.params['RasterOutputFreq']

            if nextTime >= next_intervention_time and nextTime != np.inf:
                if next_intervention_time > self.params['FinalTime']:
                    break
                self.time = next_intervention_time
                # carry out intervention update
                print("HERE", next_intervention_time)
                self.intervention_handler.update(self.all_hosts, self.time, self.all_cells)
                next_intervention_time = self.intervention_handler.next_intervention_time
            else:
                if nextTime > self.params['FinalTime']:
                    break
                # Carry out event
                self.time = nextTime
                event = self.event_handler.do_event(event_type, hostID, self.all_hosts,
                                                    self.all_cells)
                if self.params['UpdateOnAllEvents'] is True:
                    self.intervention_handler.update_on_event(event, self.all_hosts, self.time,
                                                              self.all_cells)

        self.time = self.params['FinalTime']
        end_time = time_mod.time()

        if not silent:
            print("Run {0} of {1} complete.  ".format(iteration+1, self.params['NIterations']) +
                  "Time taken: {0:.3f} seconds.".format(end_time - start_time), end="\n")

            if self.all_cells is not None:
                print("Total number of cells infected: {0}".format(np.sum([1 for x in self.all_cells if x.states["I"] > 0])))
                print("Total number of host units infected: {}".format(
                    np.sum([(x.states["C"] + x.states["I"]) for x in self.all_cells])))

        return (self.all_hosts, self.all_cells, self.run_params)

    def output_run_data(self, all_hosts, all_cells, run_params, iteration=0):
        run_data = outputdata.output_all_run_data(self, all_hosts, all_cells, run_params, iteration)
        return run_data


def run_epidemics(params, silent=False):
    all_data = []
    if params['SaveSetup']:
        run_sim = Simulator(params)
        run_sim.setup()
        run_sim.initialise(silent=silent)
    for iteration in range(params['NIterations']):
        if not params['SaveSetup']:
            run_sim = Simulator(params)
            run_sim.setup(silent=silent)
            run_sim.initialise(silent=silent)
        all_hosts, all_cells, run_params = run_sim.run_epidemic(iteration, silent=silent)
        all_data.append(
            run_sim.output_run_data(all_hosts, all_cells, run_params, iteration=iteration))

    return all_data


def main(configFile="config.ini", keyFile=False, defaultConfig=None, params_options=None,
         silent=False):
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
        config.check_params_valid(params)
        if params_options is not None:
            for key, value in params_options.items():
                # TODO Should check that these are valid additions
                params[key] = value


        params['call_params'] = copy.deepcopy(params)
        params['call_config_file'] = configFile
        params['call_module'] = str(modu)
        params['call_time'] = time_mod.strftime("%a, %d %b %Y %H:%M:%S", time_mod.localtime())
        params['call_version'] = __version__
        all_data = run_epidemics(params, silent=silent)

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
