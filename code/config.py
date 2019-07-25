"""Configuration options for the simulation.

default_config gives the sections and keys for the KEYFILE, with text descriptions."""

import os
import errno
from collections import OrderedDict
import configparser

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
        ('SimulationType', (True, "INDIVIDUAL", "Type of simulation to run.  Options are: "
                            "INDIVIDUAL, RASTER.", str)),
        ('VirtualSporulationStart', (False, None, "Distance in cells to start using VS rather than"
                                     "coupling.  Default: No VS", int)),
        ('FinalTime', (True, 10.0, "Time to stop the simulation.", float)),
        ('HostPosFile', (True, "hosts.txt", "Name of file containing host locations.  Can also "
                         "specify comma separated list of multiple files.", str)),
        ('InitCondFile', (True, "hosts_init.txt", "Name of file containing initial host states."
                          "If HostPosFile is a list, this must be a corresponding comma "
                          "separated list", str)),
        ('RegionFile', (False, None, "Name of file containing region name for each"
                        " host.  If not specified no region based data is produced.  If "
                        "HostPosFile is a list, this must be a corresponding comma separated "
                        "list", str)),
        ('KernelFile', (False, None, "Name of file containing raster kernel if SimulationType "
                        "is RASTER.", str)),
        ('NIterations', (False, 1, "Number of individual simulations to run",
                         int)),
        ('NRegions', (False, 1, "Number of distinct regions.", int)),
        ('MaxHosts', (False, 100, "Maximum number of hosts per cell in raster model.", int)),
    ])),
    ('Output', OrderedDict([
        ('OutputHostData', (False, True, "Whether to output transition times for each host.",
                            bool)),
        ('OutputEventData', (False, True, "Whether to output all individual event details.",
                             bool)),
        ('RasterOutputFreq', (True, 0.2, "How often to output raster of simulation state.  "
                              "Set to zero to suppress output.", float)),
        ('OutputFiles', (False, True, "Whether to output data to files.  "
                         "If False only python objects are returned.", bool)),
        ('OutputFileStub', (False, "output", "File path stub for output files", str)),
        ('RasterFileStub', (False, "raster_output", "File path stub for output raster files", str)),
    ])),
    ('Optimisation', OrderedDict([
        ('SaveSetup', (False, True, "Whether or not to save the initial rates and states to re-use."
                       " This can be slower for large numbers of hosts.", bool)),
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
        ('InterventionScripts', (False, None, "Comma separated list of intervention module "
                                 "scripts.  Each must include Intervention class with required "
                                 "members.", str)),
        ('InterventionUpdateFrequencies', (False, None, "Comma separated list specifying how "
                                           "often to update each intervention module.", str)),
        ('UpdateOnAllEvents', (False, False, "Whether or not to update the intervention classes "
                               "after every event.", bool)),
    ])),
])


def write_keyfile(filename="KEYFILE.ini"):
    """Generate keyfile detailing all parameter options."""

    with open(filename, "w") as outfile:
        outfile.write("# KEYFILE giving all parameter options that can be specified"
                      " in the configuration file.\n")
        for section in default_config:
            outfile.write("\n[" + section + "]\n")
            for key in default_config[section]:
                val = default_config[section][key]
                outfile.write(key + " = " + str(val[1]))
                if val[0] is False:
                    outfile.write("\n  # Optional")

                if val[2] is not None:
                    outfile.write("\n  # " + val[2] + "\n")


def write_default_config(filename="config.ini"):
    """Write config file using all default values"""

    with open(filename, "w") as outfile:
        for section in default_config:
            outfile.write("[" + section + "]\n")
            for key in default_config[section]:
                val = default_config[section][key]
                if val[0] is True:
                    outfile.write(key + " = " + str(val[1]) + "\n")
            outfile.write("\n")


def read_parser(parser):
    """Read configuration from config file parser."""

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
                    if val == "None":
                        return_dict[key] = None
                    else:
                        return_dict[key] = def_val[3](val)
            else:
                if def_val[0] is True:
                    raise configparser.NoOptionError(key, section)
                else:
                    return_dict[key] = def_val[1]

    return return_dict


def read_config_file(filename="config.ini"):
    """Read configuration file."""
    
    parser = configparser.ConfigParser()
    if os.path.exists(filename):
        parser.read(filename)
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                filename)

    return_dict = read_parser(parser)

    return return_dict


def read_config_string(config_string):
    """ Read configuration from string (for reading from log file)."""

    parser = configparser.ConfigParser()
    parser.read_string(config_string)

    return_dict = read_parser(parser)

    return return_dict


def get_parser_from_params(params):
    """Create parser for given parameters, for writing to log file."""

    parser = configparser.ConfigParser()
    parser.optionxform = str

    for section in default_config:
        parser.add_section(section)
        for key in default_config[section]:
            parser[section][key] = str(params[key])

    return parser


def check_params_valid(params):
    """Check parameter values extracted are valid."""

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
