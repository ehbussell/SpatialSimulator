"""Tools for handling output data files."""

from ..code import config
import pandas as pd


def extract_output_data(output_file_stub):
    """From an output_file_stub, extract all data as pandas dataframes."""

    params = extract_params(log_file=output_file_stub+".log")

    if params['OutputFiles'] is False:
        raise FileNotFoundError("No Output Files!")

    return_data = []
    niters = params['NIterations']
    nregions = params['NRegions']

    for run in range(niters):
        run_data = {}

        if params['OutputHostData'] is True:
            # extract host data
            filename = output_file_stub + "_hosts_" + str(run) + ".csv"
            host_data = pd.read_csv(filename)
            run_data['host_data'] = host_data

        if params['OutputEventData'] is True:
            # extract event data
            filename = output_file_stub + "_events_" + str(run) + ".csv"
            event_data = pd.read_csv(filename)
            run_data['event_data'] = event_data

        if params['SummaryOutputFreq'] > 0:
            # extract summary data
            summary_data = {}
            for region in range(nregions):
                filename = output_file_stub + "_DPC_region" + str(region) + "_" + str(run) + ".csv"
                summary_data['Region' + str(region)] = pd.read_csv(filename)
            run_data['summary_data'] = summary_data

        return_data.append(run_data)

    return return_data


def extract_params(config_file=None, log_file=None):
    """Extract the run parameters used, either from a config file or a log file."""

    if config_file is None and log_file is None:
        raise ValueError("Must specify a config_file of log_file!")

    if config_file is not None and log_file is not None:
        raise ValueError("Cannot specify both a config_file and a log_file!")

    if config_file is not None:
        # read params from config_file
        params = config.read_config_file(config_file)
    else:
        # read params from log file
        with open(log_file, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "Configuration File Used" in line:
                    config_str = "".join(lines[i+1:])
                    break

        params = config.read_config_string(config_str)

    return params
