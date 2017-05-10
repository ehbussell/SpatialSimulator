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
