"""Tools for handling output data files."""

from ..code import config
from RasterModel import raster_model
import simulator_utils
import pandas as pd
import numpy as np


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


def create_raster_run(output_file_stub, timestep=0.01, target_raster=None,
                      ignore_outside_raster=False, max_hosts=100):
    params = extract_params(log_file=output_file_stub+".log")
    all_data = extract_output_data(output_file_stub)

    if target_raster is None:
        host_raster = simulator_utils.read_raster(params['HostPosFile'].split(",")[0])
        nrows, ncols = host_raster.array.shape
    else:
        host_raster = target_raster
        nrows, ncols = host_raster.array.shape

    host_map = {}
    state = np.zeros((nrows*ncols, 2))

    s_dict = {"Cell"+str(i): [] for i in range(nrows*ncols)}
    i_dict = {"Cell"+str(i): [] for i in range(nrows*ncols)}
    f_dict = {"Cell"+str(i): [] for i in range(nrows*ncols)}
    times = []

    # Construct initial state
    for i, x in enumerate(all_data[0]['host_data'].columns):
        if "timeEnter" in x:
            host_start_idx = i
            break

    for index, host in all_data[0]['host_data'].iterrows():
        # find cell index
        cell = get_cell(host, host_raster.header_vals)
        if cell == -1:
            if ignore_outside_raster:
                continue
            else:
                raise ValueError("Host not in raster!")
        host_map[host['hostID']] = cell
        init_state = host['initial_state']
        if init_state == "S":
            state[cell, 0] += 1
        elif init_state == "I":
            state[cell, 1] += 1
        else:
            raise ValueError("Not S or I!")

    save_state(times, s_dict, i_dict, f_dict, 0.0, state, nrows*ncols)
    
    print("Initial state complete")

    next_time = timestep

    for index, event in all_data[0]['event_data'].iterrows():
        if event['time'] > next_time:
            save_state(times, s_dict, i_dict, f_dict, next_time, state, nrows*ncols)
            next_time += timestep
        if event['hostID'] in host_map:
            cell = host_map[event['hostID']]
            state[cell, 0] -= 1
            state[cell, 1] += 1
    
    save_state(times, s_dict, i_dict, f_dict, next_time, state, nrows*ncols)

    s_dict['time'] = times
    i_dict['time'] = times
    f_dict['time'] = times

    results_s = pd.DataFrame(s_dict)
    results_i = pd.DataFrame(i_dict)
    results_f = pd.DataFrame(f_dict)

    model_params = {
        'dimensions': (nrows, ncols),
        'max_hosts': max_hosts
    }

    return raster_model.RasterRun(model_params, (results_s, results_i, results_f))


def get_cell(host_row, raster_header):
    x = host_row['posX']
    y = host_row['posY']

    cellsize = raster_header['cellsize']

    col = (x - raster_header['xllcorner']) / cellsize
    row = ((raster_header['nrows'] * cellsize) -
           (y - raster_header['yllcorner'])) / cellsize

    if col < 0 or col >= raster_header['ncols']:
        return -1

    if row < 0 or row >= raster_header['nrows']:
        return -1

    cell_id = int(col) + (int(row) * raster_header['ncols'])

    return cell_id


def save_state(times, s_dict, i_dict, f_dict, time, state, ncells):
    for i in range(ncells):
        s_dict["Cell"+str(i)].append(state[i, 0])
        i_dict["Cell"+str(i)].append(state[i, 1])
        f_dict["Cell"+str(i)].append(0.0)
    times.append(time)
