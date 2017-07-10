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
        self.jump_times = {}

        if state is not None:
            self.state = state
            self.jump_times[str(state)] = 0.0

    def __repr__(self):
        repr_str = "Host(" + str(self.x) + ", " + str(self.y) + ", '"
        repr_str += str(self.state) + "', " + str(self.reg) + ")"

        return repr_str


def read_host_files(host_pos_files, init_cond_files, region_files):

    assert len(host_pos_files) == len(init_cond_files), "Number of input files do not match!"

    if region_files is not None:
        region_files = region_files.split(",")
        assert len(host_pos_files) == len(region_files), "Number of input files do not match!"
    else:
        region_files = [None for _ in host_pos_files]

    all_hosts = []

    for i, host_file in enumerate(host_pos_files):
        hosts = read_host_file(host_file, default_region=i, hostID_start=len(all_hosts))
        read_init_cond(hosts, init_cond_files[i])
        if region_files[i] is not None:
            read_regions(hosts, region_files[i])
        all_hosts.extend(hosts)

    return all_hosts


def read_host_file(filename, default_region=0, hostID_start=0):
    with open(filename, "r") as f:
        nhosts = int(f.readline())

        all_hosts = []

        for i in range(hostID_start, hostID_start + nhosts):
            x, y = f.readline().split()
            all_hosts.append(Host(float(x), float(y), hostID=i, reg=default_region))

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
