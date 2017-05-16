import simulator
import subprocess
import numpy as np
from scipy.optimize import brentq


def test_nonspatial(nruns=100):
    """
    -   Setup simulator with non-spatial kernel
    -   Run simulator N times
    -   Calculate proportion minor epidemics and assert approx equal to 1/R0
    -   Calculate limiting fraction infected from 1-z=exp(-R0z)
        and compare with mean of major epidemics
    -   Calculate stdev of major epidemic size and test percentage aggreeing with this
    -   Plot overlapping histograms
    """

    # R0 = 1.5
    #
    # p_minor = 1/R0
    #
    #
    # def size_func(z):
    #     return (1 - z - np.exp(-R0*z))
    #
    #
    # frac_inf_major = brentq(size_func, 10e-10, 1)


def test_spatial(nruns=100):
    """
    -   Setup simulator to match a preset Webidemics configuration
    -   Run simulator N times
    -   Compare distribution of final sizes using scipy.stats.anderson_ksamp
    """


def test_kernel():
    """
    -   Cache kernel for random landscape
    -   Compare with outputted kernel from Webidemics code
    """

    # Setup Webidemics
    simulator.config.gen_rand_landscape("test/randomHosts100.txt", 100)

    subprocess.run(["test/Webidemics/WebidemicsTest.exe",
                    "hostLocs=test/randomHosts100.txt",
                    "outStub=test/kernel_test_output",
                    "printKernel=1"], stdout=subprocess.DEVNULL)

    params = simulator.config.read_config_file("test/kernel_test.ini")

    test_sim = simulator.Simulator(params)
    test_sim.setup(silent=True)

    test_passed = True

    with open("test/kernel_test_output_KernelValues.txt", "r") as f:
        for line in f:
            i, j, k_ij = line.split()
            i = int(i)
            j = int(j)
            k_ij = float(k_ij)

            if k_ij != np.round(test_sim.params['kernel_vals'][i, j], 10):
                print(i, j, k_ij, np.round(test_sim.params['kernel_vals'][i, j], 10))
                test_passed = False

    return test_passed


if __name__ == "__main__":
    print("Running test_kernel...", end="")
    test1_result = test_kernel()
    if test1_result is True:
        print("PASSED")
    else:
        print("FAILED")

    print("Running test_nonspatial...", end="")
    test2_result = test_nonspatial()
    if test2_result is True:
        print("PASSED")
    else:
        print("FAILED")

    print("Running test_spatial...", end="")
    test3_result = test_spatial()
    if test3_result is True:
        print("PASSED")
    else:
        print("FAILED")
