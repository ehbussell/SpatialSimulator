"""Tools for handling output data files."""

from ..code import config
from RasterModel import raster_model
import raster_tools
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


def create_raster_runs(output_file_stub, timestep=0.01, target_raster=None,
                       ignore_outside_raster=False, max_hosts=100):
    params = extract_params(log_file=output_file_stub+".log")
    all_data = extract_output_data(output_file_stub)

    if target_raster is None:
        host_raster = raster_tools.RasterData.from_file(params['HostPosFile'].split(",")[0])
        nrows, ncols = host_raster.array.shape
    else:
        host_raster = target_raster
        nrows, ncols = host_raster.array.shape

    all_raster_runs = []

    for run_number, run in enumerate(all_data):
        host_map = {}
        state = np.zeros((nrows*ncols, 2))

        s_dict = {"Cell"+str(i): [] for i in range(nrows*ncols)}
        i_dict = {"Cell"+str(i): [] for i in range(nrows*ncols)}
        f_dict = {"Cell"+str(i): [] for i in range(nrows*ncols)}
        times = []

        # Construct initial state
        for i, x in enumerate(run['host_data'].columns):
            if "timeEnter" in x:
                host_start_idx = i
                break

        for index, host in run['host_data'].iterrows():
            # find cell index
            cell = raster_tools.find_position_in_raster((host['posX'], host['posY']),
                                                        host_raster.header_vals)
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

        print("initial state complete")

        next_time = timestep

        for index, event in run['event_data'].iterrows():
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

        results_s = pd.DataFrame(s_dict).sort_index(axis=1, inplace=True)
        results_i = pd.DataFrame(i_dict).sort_index(axis=1, inplace=True)
        results_f = pd.DataFrame(f_dict).sort_index(axis=1, inplace=True)

        model_params = {
            'dimensions': (nrows, ncols),
            'max_hosts': max_hosts
        }

        all_raster_runs.append(
            raster_model.RasterRun(model_params, (results_s, results_i, results_f)))

        print("Run {0} of {1} complete".format(run_number+1, len(all_data)))

    return all_raster_runs


def save_state(times, s_dict, i_dict, f_dict, time, state, ncells):
    for i in range(ncells):
        s_dict["Cell"+str(i)].append(state[i, 0])
        i_dict["Cell"+str(i)].append(state[i, 1])
        f_dict["Cell"+str(i)].append(0.0)
    times.append(time)


def create_cell_data(output_file_stub, target_raster=None, ignore_outside_raster=False):
    """Generate dictionary for each run with cell states and events."""

    # TODO handle target raster correctly

    print("Reading in data...")

    params = extract_params(log_file=output_file_stub+".log")
    all_data = extract_output_data(output_file_stub)

    if target_raster is None:
        host_raster = raster_tools.RasterData.from_file(params['HostPosFile'].split(",")[0])
        nrows, ncols = host_raster.array.shape
    else:
        host_raster = target_raster
        nrows, ncols = host_raster.array.shape

    all_cell_data = []

    print("Extracting cell data...")

    for run_number, run in enumerate(all_data):
        host_map = {}
        states = list(params['Model'])
        state_map = {state: i for i, state in enumerate(states)}

        cell_data = {
            i: np.array([[0.0] + [0]*len(states)]) for i in range(nrows*ncols)
        }

        values = zip(run['host_data']['hostID'].values,
                     run['host_data']['cell'].values,
                     run['host_data']['initial_state'].values)

        for hostID, cell, initState in values:
            # find cell index
            host_map[hostID] = cell
            cell_data[cell][0, state_map[initState]+1] += 1

        values = zip(run['event_data']['time'].values,
                     run['event_data']['hostID'].values,
                     run['event_data']['oldState'].values,
                     run['event_data']['newState'].values)

        for time, hostID, old, new in values:
            cell = host_map[hostID]
            new_row = np.copy(cell_data[cell][-1])
            new_row[0] = time
            new_row[state_map[old]+1] -= 1
            new_row[state_map[new]+1] += 1
            cell_data[cell] = np.vstack([cell_data[cell], new_row])

        all_cell_data.append(cell_data)
        
        print("Run {0} of {1} complete".format(run_number+1, len(all_data)))

    print("Extraction complete.")

    return all_cell_data