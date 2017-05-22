import configparser
import os
import errno
from collections import OrderedDict
from host import Host
import numpy as np

default_config = OrderedDict([
    ('Epidemiology', OrderedDict([
        ('Model', (True, "SIR", "Compartmental model to use. SECIR or SEDIR "
                   "or a subset of one of these", str)),
        ('InfRate', (True, 1.0, "Value of the infection rate.", float)),
        ('EAdvRate', (False, 1.0, "E->nextState transition rate.  Required if E in model", float)),
        ('CAdvRate', (False, 1.0, "C->nextState transition rate.  Required if C in model", float)),
        ('DAdvRate', (False, 1.0, "D->nextState transition rate.  Required if D in model", float)),
        ('IAdvRate', (False, 1.0, "I->R transition rate.  Required if R in model", float)),
        ('RAdvRate', (False, 0.0, "R->S transition rate.  Required if RS in model", float)),
        ('KernelType', (True, "EXPONENTIAL", "Functional form to use for the dispersal kernel.  "
                        "Options are: EXPONENTIAL, NONSPATIAL", str)),
        ('KernelScale', (False, 1.0, "Scale parameter for the kernel.  Required if KernelType is"
                         "EXPONENTIAL", float)),
    ])),
    ('Simulation', OrderedDict([
        ('FinalTime', (True, 10.0, "Time to stop the simulation.", float)),
        ('SummaryOutputFreq', (True, 0.2, "How often to output region summary"
                               " to output file.", float)),
        ('HostFile', (True, "hosts.txt", "Name of file containing host locations.", str)),
        ('NIterations', (False, 1, "Number of individual simulations to run",
                         int)),
        ('NRegions', (False, 1, "Number of distinct regions.", int)),
    ])),
    ('Optimisation', OrderedDict([
        ('CacheKernel', (False, False, "Whether or not to cache the full "
                         "kernel at the start of the simulation", bool)),
        ('RateStructure-Infection', (False, "ratesum",
                                     "Which rate structure to use for infection events.  Options "
                                     "are: ratesum, rateinterval, ratetree, rateCR",
                                     str)),
        ('RateStructure-Advance', (False, "ratesum",
                                   "Which rate structure to use for advance events.  Options are: "
                                   "ratesum, rateinterval, ratetree, rateCR",
                                   str)),
    ])),
])


def write_keyfile(filename="KEYFILE.txt"):
    with open(filename, "w") as f:
        f.write("# KEYFILE giving all parameter options that can be specified"
                " in the configuration file.\n# KEY* indicates that key is "
                "optional.\n")
        for section in default_config:
            f.write("\n[" + section + "]\n")
            for key in default_config[section]:
                val = default_config[section][key]
                if val[0] is False:
                    f.write(key + "* = " + str(val[1]))
                else:
                    f.write(key + " = " + str(val[1]))

                if val[2] is not None:
                    f.write("\n  # " + val[2] + "\n")


def write_default_config(filename="config.ini"):
    with open(filename, "w") as f:
        for section in default_config:
            f.write("[" + section + "]\n")
            for key in default_config[section]:
                val = default_config[section][key]
                if val[0] is True:
                    f.write(key + " = " + str(val[1]) + "\n")
            f.write("\n")


def read_config_file(filename="config.ini"):
    parser = configparser.ConfigParser()
    if os.path.exists(filename):
        parser.read(filename)
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                filename)

    return_dict = {}

    for section in default_config:
        if section not in parser:
            parser.add_section(section)

        for key in default_config[section]:
            def_val = default_config[section][key]
            if key in parser[section]:
                if def_val[3] == bool:
                    val = parser[section].getboolean(key)
                    return_dict[key] = val
                else:
                    val = parser[section][key]
                    return_dict[key] = def_val[3](val)
            else:
                if def_val[0] is True:
                    raise configparser.NoOptionError(key, section)
                else:
                    return_dict[key] = def_val[1]

    return return_dict


def read_hosts(filename="hosts.txt"):
    with open(filename, "r") as f:
        nhosts = int(f.readline())

        all_hosts = []

        for i in range(nhosts):
            x, y, state, reg = f.readline().split()
            all_hosts.append(Host(float(x), float(y), state, int(reg)))

    return all_hosts


def gen_rand_landscape(filename="hosts.txt", nhosts=100, rand_infs=0, fixed_infs=None):
    all_x = np.random.random_sample(nhosts)
    all_y = np.random.random_sample(nhosts)

    host_list = list(range(nhosts))
    infected = []

    if fixed_infs is not None:
        infected = fixed_infs
        host_list = [x for x in host_list if x not in infected]

    infected += list(np.random.choice(host_list, rand_infs, replace=False))

    with open(filename, "w") as f:
        f.write(str(nhosts) + "\n")
        for i in range(nhosts):
            if i in infected:
                f.write(str(all_x[i]) + " " + str(all_y[i]) + " I 0\n")
            else:
                f.write(str(all_x[i]) + " " + str(all_y[i]) + " S 0\n")


def verify_params(params):
    # TODO code to verify a parameter dictionary is complete
    # Option to fill in default values for missing keys with/without warning
    # Should somehow check all necessary info is present i.e. AdvRates are correctly specified
    #   even though these keys are optional
    pass
