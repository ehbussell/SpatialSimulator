def output_data_hosts(all_hosts, params, filename="output_hosts.csv"):
    states = list(params['Model'])
    nhosts = len(all_hosts)

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


def output_data_events(params, filename="output_events.csv"):
    with open(filename, "w") as f:
        f.write("time,hostID,oldState,newState\n")

        for event in params['all_events']:
            event_str = ",".join([str(x) for x in event])
            f.write(event_str + "\n")


def output_data_summary(params, filename_stub="output_DPC"):
    for region in range(params['NRegions']):
        filename = filename_stub + "_region" + str(region) + ".csv"
        states = list(params['Model'])

        with open(filename, "w") as f:
            f.write("time")
            for state in states:
                f.write(","+state)
            f.write("\n")

            for entry in params['summary_dump']:
                f.write(str(entry[0]))
                for state in states:
                    f.write(","+str(entry[1][region][state]))
                f.write("\n")
