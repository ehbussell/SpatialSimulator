"""Methods for outputting run data from the simulation."""

import inspect
import pandas as pd


def output_all_run_data(parent_sim, all_hosts, run_params, iteration=0):
    return_data = {}
    filestub = parent_sim.params['OutputFileStub']

    if parent_sim.params['OutputHostData'] is True:
        host_data = output_data_hosts(all_hosts, parent_sim.params, iteration=iteration,
                                      file_stub=filestub)
        return_data['host_data'] = host_data

    if parent_sim.params['OutputEventData'] is True:
        event_data = output_data_events(parent_sim.params, run_params, iteration=iteration,
                                        file_stub=filestub)
        return_data['event_data'] = event_data

    if parent_sim.params['SummaryOutputFreq'] != 0:
        summary_data = output_data_summary(parent_sim.params, run_params, iteration=iteration,
                                           file_stub=filestub)
        return_data['summary_data'] = summary_data

    return return_data


def output_data_hosts(all_hosts, params, file_stub="output", iteration=0):
    """Output state transition times for each host."""

    states = list(params['Model']) + ["Culled"]
    nhosts = len(all_hosts)

    filename = file_stub + "_hosts_" + str(iteration) + ".csv"

    col_names = ["hostID", "posX", "posY"] + ["timeEnter"+state for state in states]
    data_dict = {key: [] for key in col_names}

    for i in range(nhosts):
        host_id = all_hosts[i].hostID
        x = all_hosts[i].x
        y = all_hosts[i].y

        data_dict['hostID'].append(host_id)
        data_dict['posX'].append(x)
        data_dict['posY'].append(y)

        for state in states:
            jump_time = all_hosts[i].jump_times.get(state, -1.0)
            data_dict["timeEnter"+state].append(jump_time)

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
    """Ouput region based DPC summary data."""

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

    with open(params['call_config_file'], "r") as f:
        log_text += f.read()

    filename = params['OutputFileStub'] + ".log"

    with open(filename, "w") as f:
        f.write(log_text)
