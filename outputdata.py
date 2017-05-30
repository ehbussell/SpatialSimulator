def output_data_hosts(all_hosts, params, file_stub="output_hosts", iteration=0):
    states = list(params['Model']) + ["Culled"]
    nhosts = len(all_hosts)

    filename = file_stub + "_" + str(iteration) + ".csv"

    with open(filename, "w") as f:
        f.write("hostID,posX,posY")

        for state in states:
            f.write(",timeEnter"+state)
        f.write("\n")

        for i in range(nhosts):
            f.write(str(i)+","+str(all_hosts[i].x)+","+str(all_hosts[i].y))
            for state in states:
                f.write(","+str(all_hosts[i].jump_times.get(state, -1.0)))
            f.write("\n")


def output_data_events(run_params, file_stub="output_events", iteration=0):
    filename = file_stub + "_" + str(iteration) + ".csv"

    with open(filename, "w") as f:
        f.write("time,hostID,oldState,newState\n")

        for event in run_params['all_events']:
            event_str = ",".join([str(x) for x in event])
            f.write(event_str + "\n")


def output_data_summary(params, run_params, file_stub="output_DPC", iteration=0):
    for region in range(params['NRegions']):
        filename = file_stub + "_region" + str(region)
        filename += "_" + str(iteration) + ".csv"

        states = list(params['Model']) + ["Culled"]

        with open(filename, "w") as f:
            f.write("time")
            for state in states:
                f.write(","+state)
            f.write("\n")

            for entry in run_params['summary_dump']:
                f.write(str(entry[0]))
                for state in states:
                    f.write(","+str(entry[1][region][state]))
                f.write("\n")
