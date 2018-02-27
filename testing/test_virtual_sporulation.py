import os
import glob
import unittest
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.stats import poisson
import raster_tools
from IndividualSimulator import simulator

class VSRatesTests(unittest.TestCase):
    """Test VS gives correct rates"""

    @classmethod
    def setUp(cls):
        # Setup and run simulation data, single hosts on lattice with exponential kernel

        cls._data_stub = os.path.join("testing", "VS_rates_sim_output")
        cls._beta_val = 10e3
        cls._scale_val = 1
        size = (21, 21)
        kernel_size = (43, 43)
        kernel_centre = [int(x/2) for x in kernel_size]

        # Create host file
        host_raster = raster_tools.RasterData(size, array=np.ones(size))
        host_file = os.path.join("testing", "VS_rates_host_test_case.txt")
        host_raster.to_file(host_file)

        # Create initial conditions files
        init_stub = os.path.join("testing", "VS_rates_init_test_case")
        host_raster.array[int(size[0]/2), int(size[1]/2)] = 0
        host_raster.to_file(init_stub + "_E.txt")
        host_raster.array = np.zeros(size)
        host_raster.to_file(init_stub + "_S.txt")
        host_raster.array[int(size[0]/2), int(size[1]/2)] = 1
        host_raster.to_file(init_stub + "_I.txt")

        # Create kernel file
        kernel_raster = raster_tools.RasterData(kernel_size, array=np.full(kernel_size, 1.0))
        kernel_raster.array[kernel_centre[0], kernel_centre[1]] = 0
        kernel_raster.to_file(init_stub + "_kernel.txt")

        cls._kernel = kernel_raster.array
        cls._size = size

        # Setup config file
        config_filename = os.path.join("testing", "VS_rates_config.ini")
        config_str = "\n[Epidemiology]\n"
        config_str += "Model = SEI\nInfRate = " + str(cls._beta_val) 
        config_str += "\nEAdvRate = 0.0\nIAdvRate = 0.0\nKernelType = RASTER\n"
        config_str += "\n[Simulation]\n"
        config_str += "SimulationType = RASTER\nFinalTime = 1\nNIterations = 1\n"
        config_str += "HostPosFile = " + host_file + "\nInitCondFile = " + init_stub + "\n"
        config_str += "KernelFile = " + init_stub + "_kernel.txt" + "\n"
        config_str += "VirtualSporulationStart = 1"
        config_str += "\n[Output]\n"
        config_str += "SummaryOutputFreq = 0\nOutputFileStub = " + cls._data_stub
        config_str += "\n[Optimisation]\n"
        config_str += "SaveSetup = False"
        with open(config_filename, "w") as outfile:
            outfile.write(config_str)

        cls._simulator = simulator.Simulator(config_file=config_filename)
        cls._simulator.setup()

    def test_VS_rates(self):

        self._simulator.initialise()

        self.assertTrue(np.allclose(self._kernel, self._simulator.params['kernel']))

        print(self._simulator.params['spore_rate'])

        count = Counter()
        n_spore_events = 0

        while True:
            totRate, event_type, hostID = self._simulator.rate_handler.get_next_event()
            if event_type != "Sporulation":
                raise ValueError("Not a sporulation event!")
            nextTime = self._simulator.time + (-1.0/totRate)*np.log(np.random.random_sample())
            if nextTime > self._simulator.params['FinalTime']:
                break
            n_spore_events += 1
            # Carry out event
            self._simulator.time = nextTime
            cell_id = self._simulator.event_handler.do_event(
                event_type, hostID, self._simulator.all_hosts, self._simulator.all_cells)
            if cell_id is not None:
                count[cell_id] += 1

        centre_pos = [int(self._size[i]/2) for i in range(2)]
        kernel_shape = self._kernel.shape

        all_n_events = []

        for key, val in count.items():

            host_pos = np.unravel_index(key, self._size)
            rel_pos = (host_pos[0] - centre_pos[0], host_pos[1] - centre_pos[1])
            kernel_pos = [rel_pos[i] + int(kernel_shape[i]/2) for i in range(2)]
            expected = self._beta_val*self._kernel[kernel_pos[0], kernel_pos[1]]
            all_n_events.append(val)

        mu = expected
        x = np.arange(poisson.ppf(0.001, mu), poisson.ppf(0.999, mu))

        plt.style.use("ggplot")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(all_n_events, normed=True, bins=30, color="blue", alpha=0.3)
        ax.plot(x, poisson.pmf(x, mu), 'r--', lw=2, label='poisson pmf')
        ax.set_xlabel("Number of Spores Landing")
        ax.set_ylabel("Frequency")
        ax.set_title("Virtual Sporulation Test Results")
        fig.savefig(os.path.join("testing", "NEventsHist.png"))

        bins = list(x) + [x[-1]+0.5]
        f_obs = np.histogram(all_n_events, bins, normed=True)[0]
        f_exp = poisson.pmf(x, mu)

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)


    @classmethod
    def tearDownClass(cls):
        os.remove(os.path.join("testing", "VS_rates_config.ini"))
        os.remove(os.path.join("testing", "VS_rates_host_test_case.txt"))
        for file in glob.glob(os.path.join("testing", "VS_rates_init_test_case_*")):
            os.remove(file)
