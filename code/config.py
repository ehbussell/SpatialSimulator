"""Configuration options for the simulation.

default_config gives the sections and keys for the KEYFILE, with text descriptions."""

import configparser
import os
import errno
from collections import OrderedDict
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
        ('HostPosFile', (True, "hosts.txt", "Name of file containing host locations.", str)),
        ('InitCondFile', (True, "hosts_init.txt", "Name of file containing initial host states.",
                          str)),
        ('RegionFile', (False, None, "Name of file containing region name for each"
                        " host.  If not specified no region based data is produced", str)),
        ('NIterations', (False, 1, "Number of individual simulations to run",
                         int)),
        ('NRegions', (False, 1, "Number of distinct regions.", int)),
    ])),
    ('Output', OrderedDict([
        ('OutputHostData', (False, True, "Whether to output transition times for each host.",
                            bool)),
        ('OutputEventData', (False, True, "Whether to output all individual event details.",
                             bool)),
        ('SummaryOutputFreq', (True, 0.2, "How often to output region summary.  "
                               "Set to zero to suppress output.", float)),
        ('OutputFiles', (False, True, "Whether to output data to files.  "
                         "If False only python objects are returned.", bool)),
        ('OutputFileStub', (False, "output", "File path stub for output files", str)),
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
    ('Interventions', OrderedDict([
        ('InterventionScript', (False, None, "Intervention module script.  Must include "
                                "create_interventions function that returns an intervention class",
                                str)),
        ('InterventionUpdateFrequency', (False, 0.0, "How often to carry out intervention update.",
                                         float)),
        ('UpdateOnAllEvents', (False, False, "Whether or not to update the intervention class "
                               "after every event.", bool)),
        ('ContinuousRemovalRate', (False, 0.0, "Rate at which infected hosts are culled.",
                                   float)),
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


def read_parser(parser):
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


def read_config_file(filename="config.ini"):
    parser = configparser.ConfigParser()
    if os.path.exists(filename):
        parser.read(filename)
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                filename)

    return_dict = read_parser(parser)

    return return_dict


def read_config_string(config_string):
    parser = configparser.ConfigParser()
    parser.read_string(config_string)

    return_dict = read_parser(parser)

    return return_dict


def check_params_valid(params):
    # TODO Should somehow check all necessary info is present i.e. AdvRates are correctly specified
    #   even though these keys are optional
    for section in default_config:
        for key, def_val in default_config[section].items():
            if key in params:
                if isinstance(params[key], def_val[3]) is False and params[key] is not None:
                    raise TypeError("Incorrect type for key {0}.  Expected: {1}, Found: "
                                    "{2}".format(key, def_val[3], type(params[key])))
            else:
                raise KeyError("Missing parameter key {0}".format(key))
