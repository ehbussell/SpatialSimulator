"""Methods for outputting run data from the simulation."""

import io
import numpy as np
import os
import pandas as pd
from . import config
import raster_tools

def output_all_run_data(parent_sim, all_hosts, all_cells, run_params, iteration=0):
    return_data = {}
    filestub = parent_sim.params['OutputFileStub']

    if parent_sim.params['OutputFiles'] is True:
        output_path = os.path.split(filestub)[0]
        if output_path != "":
            os.makedirs(output_path, exist_ok=True)

    if parent_sim.params['OutputHostData'] is True:
        host_data = output_data_hosts(all_hosts, all_cells, parent_sim.params, iteration=iteration,
                                      file_stub=filestub)
        return_data['host_data'] = host_data

    if parent_sim.params['OutputEventData'] is True:
        event_data = output_data_events(parent_sim.params, run_params, iteration=iteration,
                                        file_stub=filestub)
        return_data['event_data'] = event_data

    # if parent_sim.params['SummaryOutputFreq'] != 0:
    #     raise NotImplementedError
    #     summary_data = output_data_summary(parent_sim.params, run_params, iteration=iteration,
    #                                        file_stub=filestub)
    #     return_data['summary_data'] = summary_data

    if parent_sim.params['InterventionScripts'] is not None:
        control_data = parent_sim.intervention_handler.output(iteration=iteration,
                                                              file_stub=filestub)

        return_data['control_data'] = control_data

    return return_data


def output_raster_data(parent_sim, time=None, iteration=None, states=None):
    """Output current state raster"""

    all_cells = parent_sim.all_cells
    header = parent_sim.params['header']

    if states is None:
        states = list(parent_sim.params['Model']) + ["Culled"]

    output_stub = parent_sim.params['RasterFileStub']

    if iteration is not None:
        iterstub = "_" + str(iteration)
    else:
        iterstub = ""

    if time is not None:
        timestub = "_" + str(time)
    else:
        timestub = ""

    for state in states:
        cell_state = np.full((header['nrows'], header['ncols']), header['NODATA_value'])
        for row in range(header['nrows']):
            for col in range(header['ncols']):
                cell_id = parent_sim.params['cell_map'].get((row, col), None)
                if cell_id is not None:
                    cell_state[row, col] = all_cells[cell_id].states[state]
        cell_state = cell_state.reshape((header['nrows'], header['ncols']))

        raster = raster_tools.RasterData(
            shape=(header['nrows'], header['ncols']),
            llcorner=(header['xllcorner'], header['yllcorner']),
            cellsize=header['cellsize'],
            NODATA_value=header['NODATA_value'],
            array=cell_state
        )

        raster.to_file(output_stub + iterstub + "_" + state + timestub + ".txt")

def output_data_hosts(all_hosts, all_cells, params, file_stub="output", iteration=0):
    """Output state transition times for each host."""

    states = list(params['Model']) + ["Culled"]
    nhosts = len(all_hosts)

    filename = file_stub + "_hosts_" + str(iteration) + ".csv"

    if params['SimulationType'] == "INDIVIDUAL":
        col_names = ["hostID", "posX", "posY", "region", "initial_state"] + ["timeEnter"+state for state in states]
    else:
        col_names = ["hostID", "posX", "posY", "cell", "cell_pos", "region", "initial_state"] + [
            "timeEnter" + state for state in states]

    data_dict = {key: [] for key in col_names}

    for i in range(nhosts):
        host_id = all_hosts[i].host_id
        x = all_hosts[i].xpos
        y = all_hosts[i].ypos
        region = all_hosts[i].reg
        start_state = all_hosts[i].init_state

        data_dict['hostID'].append(host_id)
        data_dict['posX'].append(x)
        data_dict['posY'].append(y)
        data_dict['region'].append(region)
        data_dict['initial_state'].append(start_state)

        if params['SimulationType'] == "RASTER":
            data_dict['cell'].append(all_hosts[i].cell_id)
            data_dict['cell_pos'].append(all_cells[all_hosts[i].cell_id].cell_position)

        for state in states:
            jump_times = []
            for jump in all_hosts[i].trans_times:
                if jump[2] == state:
                    jump_times.append(jump[0])
            data_dict["timeEnter"+state].append(jump_times)

    data_frame = pd.DataFrame(data_dict, columns=col_names)

    if params['OutputFiles'] is True:
        data_frame.to_csv(filename, index=False)

    return data_frame


def output_data_events(params, run_params, file_stub="output", iteration=0):
    """Output details on each state change event, in chronological order."""

    filename = file_stub + "_events_" + str(iteration) + ".csv"

    col_names = ["time", "hostID", "oldState", "newState"]
    data = [[] for _ in col_names]

    for event in run_params['all_events']:
        for i, x in enumerate(event):
            data[i].append(x)

    data_dict = {col_names[i]: data[i] for i in range(len(col_names))}
    data_frame = pd.DataFrame(data_dict, columns=col_names)

    if params['OutputFiles'] is True:
        data_frame.to_csv(filename, index=False)

    return data_frame


def output_data_summary(params, run_params, file_stub="output", iteration=0):
    """Output region based DPC summary data."""

    # TODO re-implement now summary dump removed - construct from host data?

    data_dict = {}
    for region in range(params['NRegions']):
        filename = file_stub + "_DPC_region" + str(region) + "_" + str(iteration) + ".csv"

        states = list(params['Model']) + ["Culled"]
        col_names = ["time"] + states
        region_data_dict = {key: [] for key in col_names}

        for entry in run_params['summary_dump']:
            time = entry[0]
            region_data_dict['time'].append(time)
            for state in states:
                state_num = entry[1][region][state]
                region_data_dict[state].append(state_num)

        data_dict['Region'+str(region)] = pd.DataFrame(region_data_dict, columns=col_names)

        if params['OutputFiles'] is True:
            data_dict['Region'+str(region)].to_csv(filename, index=False)

    return data_dict


def output_log_file(params):
    log_text = "Individual Simulator Log File\n" + "#"*35 + "\n\n"
    log_text += "Simulator called from: " + params['call_module'] + "\n\n"
    log_text += "Version: " + params['call_version'] + "\n\n"
    log_text += "Run: " + params['call_time'] + "\n\n"
    log_text += "Configuration File Used\n" + "#"*35 + "\n\n"

    parser = config.get_parser_from_params(params['call_params'])
    config_txt = io.StringIO()
    parser.write(config_txt)

    log_text += config_txt.getvalue()

    filename = params['OutputFileStub'] + ".log"

    with open(filename, "w") as outfile:
        outfile.write(log_text)
