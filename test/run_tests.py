import simulator
import subprocess
import glob
import os
import argparse
import numpy as np
from scipy.optimize import brentq
from scipy.stats import anderson_ksamp
import matplotlib.pyplot as plt

plt.style.use("ggplot")
graphParams = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 8,
   'xtick.labelsize': 8,
   'ytick.labelsize': 8,
   'text.usetex': False,
   'figure.figsize': [6, 4]
}


def test_nonspatial(nruns=10000, otherOptions=None):
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

    if otherOptions is not None:
        addOtherOptions(params, otherOptions)

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

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.hist(final_rs, bins=50, normed=True)
    ax1.set_xlabel("Final Number Removed")
    ax1.set_ylabel("Probability")
    fig.savefig("test/nonspatial_test_finalSizeDistrib", dpi=300)

    return test_passed


def test_spatial(nruns=1000, otherOptions=None):
    """
    -   Setup simulator to match a preset Webidemics configuration
    -   Run simulator N times
    -   Compare distribution of final sizes using scipy.stats.anderson_ksamp
    """

    # Setup Webidemics
    simulator.config.gen_rand_landscape("test/randomHosts1000.txt", 1000, fixed_infs=[0, 1])

    subprocess.run(["test/Webidemics/WebidemicsTest.exe",
                    "hostLocs=test/randomHosts1000.txt",
                    "outStub=test/spatial_test_output",
                    "numIts="+str(nruns)], stdout=subprocess.DEVNULL)

    # Setup simulator
    params = simulator.config.read_config_file("test/spatial_test.ini")
    params['NIterations'] = nruns

    if otherOptions is not None:
        addOtherOptions(params, otherOptions)

    test_sim = simulator.Simulator(params)
    test_sim.setup(silent=True)

    Webidemics_final_rs = []
    test_sim_final_rs = []

    for i in range(nruns):
        all_hosts, run_params = test_sim.run_epidemic(i)
        test_sim_final_rs.append(run_params['summary_dump'][-1][1][0]['R'])

        with open("test/spatial_test_output_end_"+str(i)+".txt", "r") as f:
            line = f.readline()
            Webidemics_final_rs.append(int(line.split()[5]))

    bins = np.linspace(0, 1000, 100)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.hist(test_sim_final_rs, bins, alpha=0.5, label="Simulator", normed=True)
    ax1.hist(Webidemics_final_rs, bins, alpha=0.5, label="Webidemics", normed=True)
    ax1.legend(loc="upper right")
    ax1.set_xlabel("Final Number Removed")
    ax1.set_ylabel("Probability")
    fig.savefig("test/spatial_test_finalSizeDistrib", dpi=300)

    test = anderson_ksamp([test_sim_final_rs, Webidemics_final_rs])

    print(test)

    for f in glob.glob("test/spatial_test_output*"):
        os.remove(f)

    if test[0] > test[1][1]:
        return False
    else:
        return True


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


def test_continuous_removal(nruns=1000, otherOptions=None):
    """
    -   Setup simulator to match a preset Webidemics configuration
    -   Match I->R rate such that half comes from continuous removal intervention
    -   Run simulator N times
    -   Compare distribution of final sizes using scipy.stats.anderson_ksamp
    """

    # Setup Webidemics
    simulator.config.gen_rand_landscape("test/randomHosts1000.txt", 1000, fixed_infs=[0, 1])

    subprocess.run(["test/Webidemics/WebidemicsTest.exe",
                    "hostLocs=test/randomHosts1000.txt",
                    "outStub=test/cont_removal_test_output",
                    "numIts="+str(nruns)], stdout=subprocess.DEVNULL)

    # Setup simulator
    params = simulator.config.read_config_file("test/cont_removal_test.ini")
    params['NIterations'] = nruns

    if otherOptions is not None:
        addOtherOptions(params, otherOptions)

    test_sim = simulator.Simulator(params)
    test_sim.setup(silent=True)

    Webidemics_final_rs = []
    test_sim_final_rs = []

    for i in range(nruns):
        all_hosts, run_params = test_sim.run_epidemic(i)
        test_sim_final_rs.append(run_params['summary_dump'][-1][1][0]['R'] +
                                 run_params['summary_dump'][-1][1][0]['Culled'])

        with open("test/cont_removal_test_output_end_"+str(i)+".txt", "r") as f:
            line = f.readline()
            Webidemics_final_rs.append(int(line.split()[5]))

    bins = np.linspace(0, 1000, 100)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.hist(test_sim_final_rs, bins, alpha=0.5, label="Simulator", normed=True)
    ax1.hist(Webidemics_final_rs, bins, alpha=0.5, label="Webidemics", normed=True)
    ax1.legend(loc="upper right")
    ax1.set_xlabel("Final Number Removed")
    ax1.set_ylabel("Probability")
    fig.savefig("test/cont_removal_test_finalSizeDistrib", dpi=300)

    test = anderson_ksamp([test_sim_final_rs, Webidemics_final_rs])

    print(test)

    for f in glob.glob("test/cont_removal_test_output*"):
        os.remove(f)

    if test[0] > test[1][1]:
        return False
    else:
        return True


def addOtherOptions(params, otherOptions):
    for pair in otherOptions:
        key, val = pair.split("=")
        type_val = None
        for section in simulator.config.default_config:
            if key in simulator.config.default_config[section]:
                type_val = simulator.config.default_config[section][key][3]

        if type_val is not None:
            val = type_val(val)
        else:
            raise ValueError("Unknown key=val pair!")

        print(key, val)

        params[key] = val


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run tests on simulator.")
    parser.add_argument("Tests", help="Test numbers to run.  " +
                        "If blank all tests are carried out", nargs="*",
                        default=[0], type=int)
    parser.add_argument("-s", "--Spatial_nRuns",
                        help="Number of runs to carry out in spatial test.", default=1000,
                        type=int)
    parser.add_argument("-n", "--NonSpatial_nRuns",
                        help="Number of runs to carry out in non-spatial test.", default=10000,
                        type=int)
    parser.add_argument("-c", "--Cont_nRuns",
                        help="Number of runs to carry out in continuous removal test.",
                        default=1000, type=int)
    parser.add_argument("-o", "--otherOptions", help="Further key=val pairs to pass to params.",
                        nargs="*", default=None, type=str)
    args = parser.parse_args()

    test_list = args.Tests

    if 1 in test_list or 0 in test_list:
        print("Running test_kernel...")
        test1_result = test_kernel()
        if test1_result is True:
            print("...test_kernel PASSED")
        else:
            print("...test_kernel FAILED")

    if 2 in test_list or 0 in test_list:
        print("Running test_nonspatial...")
        test2_result = test_nonspatial(args.NonSpatial_nRuns, args.otherOptions)
        if test2_result is True:
            print("...test_nonspatial PASSED")
        else:
            print("...test_nonspatial FAILED")

    if 3 in test_list or 0 in test_list:
        print("Running test_spatial...")
        test3_result = test_spatial(args.Spatial_nRuns, args.otherOptions)
        if test3_result is True:
            print("...test_spatial PASSED")
        else:
            print("...test_spatial FAILED")

    if 4 in test_list or 0 in test_list:
        print("Running test_continuous_removal...")
        test4_result = test_continuous_removal(args.Cont_nRuns, args.otherOptions)
        if test4_result is True:
            print("...test_continuous_removal PASSED")
        else:
            print("...test_continuous_removal FAILED")
