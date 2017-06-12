"""Methods for handling host detail initialisation, including all reading of host files."""


class Host(object):
    """All stored data for an individual host.

    Attributes:
        x:          X position of host
        y:          Y position of host
        state:      Current state of host (S, I, R etc)
        reg:        Region containing host
        hostID:     ID code unique to this host
        jump_times: Dictionary of entry times into each state compartment that has been visited.
    """

    def __init__(self, x, y, state=None, reg=0, hostID=None):
        self.x = x
        self.y = y
        self.reg = reg
        self.hostID = hostID

        if state is not None:
            self.state = state
            self.jump_times = {str(state): 0.0}

    def __repr__(self):
        repr_str = "Host(" + str(self.x) + ", " + str(self.y) + ", '"
        repr_str += str(self.state) + "', " + str(self.reg) + ")"

        return repr_str


def read_host_files(host_pos_file, init_cond_file, region_file):
    all_hosts = read_host_file(host_pos_file)

    read_init_cond(all_hosts, init_cond_file)
    if region_file is not None:
        read_regions(all_hosts, region_file)

    return all_hosts


def read_host_file(filename):
    with open(filename, "r") as f:
        nhosts = int(f.readline())

        all_hosts = []

        for i in range(nhosts):
            x, y = f.readline().split()
            all_hosts.append(Host(float(x), float(y), i))

    return all_hosts


def read_init_cond(all_hosts, filename):
    with open(filename, "r") as f:
        nhosts = int(f.readline())

        for i in range(nhosts):
            state = f.readline().strip()
            all_hosts[i].state = state
            all_hosts[i].jump_times[str(state)] = 0.0


def read_regions(all_hosts, filename):
    with open(filename, "r") as f:
        nhosts = int(f.readline())

        for i in range(nhosts):
            region = int(f.readline().strip())
            all_hosts[i].reg = region
