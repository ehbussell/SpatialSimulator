import pandas as pd


def output_all_run_data(parent_sim, all_hosts, run_params, iteration=0):
    return_data = {}

    if parent_sim.params['OutputHostData'] is True:
        host_data = output_data_hosts(all_hosts, parent_sim.params, iteration=iteration)
        return_data['host_data'] = host_data

    if parent_sim.params['OutputEventData'] is True:
        event_data = output_data_events(run_params, iteration=iteration)
        return_data['event_data'] = event_data

    if parent_sim.params['SummaryOutputFreq'] != 0:
        summary_data = output_data_summary(parent_sim.params, run_params, iteration=iteration)
        return_data['summary_data'] = summary_data

    return return_data


def output_data_hosts(all_hosts, params, file_stub="output_hosts", iteration=0):
    states = list(params['Model']) + ["Culled"]
    nhosts = len(all_hosts)

    filename = file_stub + "_" + str(iteration) + ".csv"

    col_names = ["hostID", "posX", "posY"] + ["timeEnter"+state for state in states]
    data_dict = {key: [] for key in col_names}

    with open(filename, "w") as f:
        f.write(",".join(col_names))
        f.write("\n")

        for i in range(nhosts):
            host_id = all_hosts[i].hostID
            x = all_hosts[i].x
            y = all_hosts[i].y

            data_dict['hostID'].append(host_id)
            data_dict['posX'].append(x)
            data_dict['posY'].append(y)

            f.write(str(host_id)+","+str(x)+","+str(y))

            for state in states:
                jump_time = all_hosts[i].jump_times.get(state, -1.0)
                data_dict["timeEnter"+state].append(jump_time)
                f.write(","+str(jump_time))

            f.write("\n")

    data_frame = pd.DataFrame(data_dict, columns=col_names)

    return data_frame


def output_data_events(run_params, file_stub="output_events", iteration=0):
    filename = file_stub + "_" + str(iteration) + ".csv"

    col_names = ["time", "hostID", "oldState", "newState"]
    data = [[] for _ in col_names]

    with open(filename, "w") as f:
        f.write(",".join(col_names))

        for event in run_params['all_events']:
            event_str = ",".join([str(x) for x in event])
            for i, x in enumerate(event):
                data[i].append(x)
            f.write(event_str + "\n")

    data_dict = {col_names[i]: data[i] for i in range(len(col_names))}
    data_frame = pd.DataFrame(data_dict, columns=col_names)

    return data_frame


def output_data_summary(params, run_params, file_stub="output_DPC", iteration=0):
    data_dict = {}
    for region in range(params['NRegions']):
        filename = file_stub + "_region" + str(region)
        filename += "_" + str(iteration) + ".csv"

        states = list(params['Model']) + ["Culled"]
        col_names = ["time"] + states
        region_data_dict = {key: [] for key in col_names}

        with open(filename, "w") as f:
            f.write("time")
            for state in states:
                f.write(","+state)
            f.write("\n")

            for entry in run_params['summary_dump']:
                time = entry[0]
                f.write(str(time))
                region_data_dict['time'].append(time)
                for state in states:
                    state_num = entry[1][region][state]
                    f.write(","+str(state_num))
                    region_data_dict[state].append(state_num)
                f.write("\n")

        data_dict['Region'+str(region)] = pd.DataFrame(region_data_dict, columns=col_names)

    return data_dict
