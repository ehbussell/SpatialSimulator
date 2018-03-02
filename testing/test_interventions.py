import unittest
import os
import pdb
import glob
from scipy.stats import expon, kstest
import numpy as np
import matplotlib.pyplot as plt
import raster_tools
import IndividualSimulator

full_lambda = 0.2

class ContinuousIntervention():
    """Test intervention implementing continuous removal of infected hosts."""

    def __init__(self, update_freq, all_hosts, all_cells=None):
        self.update_freq = update_freq
        self.type = "CONTINUOUS"
        self.rate_size = 1
        self.rate_factor = full_lambda / 2

    def action(self, all_hosts, time, event_id, all_cells=None):
        for host in all_hosts:
            if host.state == "I":
                return [(host.host_id, "CULL")]

    def update(self, all_hosts, time, all_cells=None, after_event=None, get_rate_fn=None,
               initial=False):

        if initial:
            total_inf = 0
            for host in all_hosts:
                if host.state == "I":
                    total_inf += 1
            return [(0, total_inf)]

        if after_event is None:
            return []
        else:
            host_id, cell_id, old_state, new_state = after_event
            if old_state == "I":
                rate_change = -1.0
            if new_state == "I":
                rate_change = +1.0

            new_rate = get_rate_fn(0) + rate_change
            return [(0, new_rate)]

    def finalise(self):
        pass

    # Function to get log string to include in log file for set of simulations. Called after all
    # simulations have been carried out and finalise has been called for each.
    def output(self):
        log_str = ""
        return log_str


class TestContinuousInterventions(unittest.TestCase):
    """Test implementation of continuous interventions."""

    def setUp(self):
        # Setup and run simulation data, single infected hosts on lattice with non-spatial kernel.

        self._full_lambda = full_lambda
        self._data_stub = os.path.join("testing", "cont_intervention_sim_output")
        size = (100, 100)
        kernel_size = (3, 3)
        kernel_centre = [int(x/2) for x in kernel_size]

        # Create host file
        host_raster = raster_tools.RasterData(size, array=np.ones(size))
        host_file = os.path.join("testing", "cont_intervention_host_test_case.txt")
        host_raster.to_file(host_file)

        # Create initial conditions files
        init_stub = os.path.join("testing", "cont_intervention_init_test_case")
        host_raster.to_file(init_stub + "_I.txt")
        host_raster.array = np.zeros(size)
        host_raster.to_file(init_stub + "_S.txt")
        host_raster.to_file(init_stub + "_R.txt")

        # Create kernel file
        kernel_raster = raster_tools.RasterData(kernel_size, array=np.full(kernel_size, 1.0))
        kernel_raster.array[kernel_centre[0], kernel_centre[1]] = 0
        kernel_raster.to_file(init_stub + "_kernel.txt")

        self._kernel = kernel_raster.array
        self._size = size

        # Setup config file
        self.config_filename = os.path.join("testing", "cont_intervention_config.ini")
        config_str = "\n[Epidemiology]\n"
        config_str += "Model = SIR\nInfRate = 0"
        config_str += "\nIAdvRate = " + str(self._full_lambda/2) + "\nKernelType = RASTER\n"
        config_str += "\n[Simulation]\n"
        config_str += "SimulationType = RASTER\nFinalTime = 1000\nNIterations = 1\n"
        config_str += "HostPosFile = " + host_file + "\nInitCondFile = " + init_stub + "\n"
        config_str += "KernelFile = " + init_stub + "_kernel.txt" + "\n"
        config_str += "VirtualSporulationStart = None"
        config_str += "\n[Output]\n"
        config_str += "SummaryOutputFreq = 0\nOutputFileStub = " + self._data_stub
        config_str += "\n[Optimisation]\n"
        config_str += "SaveSetup = False"
        config_str += "\n[Interventions]\n"
        config_str += "InterventionUpdateFrequencies = None\n"
        config_str += "UpdateOnAllEvents = True\n"

        with open(self.config_filename, "w") as outfile:
            outfile.write(config_str)

    def test_cont_int_nonspatial(self):
        """Test continuous removal with non-spatial epidemic."""

        params = IndividualSimulator.code.config.read_config_file(filename=self.config_filename)
        params['InterventionScripts'] = [ContinuousIntervention]

        simulator = IndividualSimulator.Simulator(params=params)
        simulator.setup()
        simulator.initialise()
        all_hosts, all_cells, run_params = simulator.run_epidemic()

        waiting_times = [host.trans_times[1][0] for host in all_hosts]

        scale = 1.0 / self._full_lambda
        x = np.linspace(expon.ppf(0.001, scale=scale), expon.ppf(0.999, scale=scale), num=100)

        plt.style.use("ggplot")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(waiting_times, normed=True, bins=100, color="blue", alpha=0.3)
        ax.plot(x, expon.pdf(x, scale=scale), 'r--', lw=2, label='Exponential pdf')
        ax.set_xlabel("Waiting Time")
        ax.set_ylabel("Frequency")
        ax.set_title("Continuous Intervention Test Results")
        fig.savefig(os.path.join("testing", "ContInterventionHist.png"))

        ks_stat, pval = kstest(waiting_times, expon.cdf, args=(0, scale))
        self.assertGreater(
            pval, 0.1, msg="Waiting Time distribution significantly different from Exponential.")

    def tearDown(self):
        os.remove(os.path.join("testing", "cont_intervention_config.ini"))
        os.remove(os.path.join("testing", "cont_intervention_host_test_case.txt"))
        for file in glob.glob(os.path.join("testing", "cont_intervention_init_test_case_*")):
            os.remove(file)
