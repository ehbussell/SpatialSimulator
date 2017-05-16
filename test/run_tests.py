import simulator
import subprocess
import glob
import os
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

plt.style.use("ggplot")


def test_nonspatial(nruns=10000):
    """
    -   Setup simulator with non-spatial kernel
    -   Run simulator N times
    -   Calculate proportion minor epidemics and assert approx equal to 1/R0
    -   Calculate limiting fraction infected from 1-z=exp(-R0z)
        and compare with mean of major epidemics
    -   Calculate stdev of major epidemic size and test percentage aggreeing with this
    -   Plot overlapping histograms
    """

    simulator.config.gen_rand_landscape("test/randomHosts1000.txt", 1000, rand_infs=1)
    params = simulator.config.read_config_file("test/nonspatial_test.ini")
    params['NIterations'] = nruns

    final_rs = []

    test_sim = simulator.Simulator(params)
    test_sim.setup(silent=True)
    for iteration in range(params['NIterations']):
        all_hosts, run_params = test_sim.run_epidemic(iteration)
        final_rs.append(run_params['summary_dump'][-1][1][0]['R'])

    R0 = params['InfRate']*1000/params['IAdvRate']

    exp_p_minor = 1/R0

    def size_func(z):
        return (1 - z - np.exp(-R0*z))

    exp_frac_inf_major = brentq(size_func, 10e-10, 1)

    major_rs = [x for x in final_rs if x > 100]
    p_minor = (nruns - len(major_rs))/nruns
    frac_inf_major = np.sum(major_rs)/(len(major_rs)*1000)

    print(exp_p_minor, p_minor)
    print(exp_frac_inf_major, frac_inf_major)

    test_passed = True
    if np.round(p_minor, 3) != np.round(exp_p_minor, 3):
        test_passed = False
    if np.round(frac_inf_major, 3) != np.round(exp_frac_inf_major, 3):
        test_passed = False

    plt.hist(final_rs, bins=50)
    plt.show()

    return test_passed


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

    # Setup simulator
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

    for f in glob.glob("test/kernel_test_output*"):
        os.remove(f)

    return test_passed


if __name__ == "__main__":
    print("Running test_kernel...")
    test1_result = test_kernel()
    if test1_result is True:
        print("...test_kernel PASSED")
    else:
        print("...test_kernel FAILED")

    print("Running test_nonspatial...")
    test2_result = test_nonspatial()
    if test2_result is True:
        print("...test_nonspatial PASSED")
    else:
        print("...test_nonspatial FAILED")

    print("Running test_spatial...")
    test3_result = test_spatial()
    if test3_result is True:
        print("...test_spatial PASSED")
    else:
        print("...test_spatial FAILED")
