"""Methods for handling host detail initialisation, including all reading of host files."""

import numpy as np
import raster_tools


class Host(object):
    """All stored data for an individual host.

    Attributes:
        xpos:           X position of host
        ypos:           Y position of host
        state:          Current state of host (S, I, R etc)
        init_state:     Start state of host (at time t=0)
        reg:            Region containing host
        host_id:        ID code unique to this host
        trans_times:    List of state transitions.  Each entry is tuple (time, old state, new state)
    """

    def __init__(self, x, y, state=None, reg=0, host_id=None, cell_id=None):
        self.xpos = x
        self.ypos = y
        self.reg = reg
        self.host_id = host_id
        self.cell_id = cell_id
        self.trans_times = []

        if state is not None:
            self.state = state
            self.init_state = state
            self.trans_times.append((0.0, None, state))

    def __repr__(self):
        repr_str = "Host(" + str(self.xpos) + ", " + str(self.ypos) + ", '"
        repr_str += str(self.state) + "', " + str(self.reg) + str(self.host_id) + ")"

        return repr_str

    def initialise_state(self, start_state):
        """Update state at start time."""

        self.trans_times.append((0.0, None, start_state))
        self.state = start_state
        self.init_state = start_state

    def update_state(self, new_state, time):
        """Update the current state of the host, storing transition time."""

        self.trans_times.append((time, self.state, new_state))
        self.state = new_state


class Cell(object):
    """Raster cell object holding multiple hosts."""

    def __init__(self, cell_position, hosts=None, cell_id=None):
        self.cell_id = cell_id
        self.cell_position = cell_position
        all_states = ["S", "E", "C", "D", "I", "R", "Culled"]
        self.states = {state: 0 for state in all_states}
        self.susceptibility = 1
        self.infectiousness = 1

        if hosts is not None:
            self.hosts = hosts
            for host in self.hosts:
                if host.state is not None:
                    self.states[host.state] = self.states[host.state] + 1
        else:
            self.hosts = []

    def __repr__(self):
        repr_str = "Cell(" + str(self.cell_position) + ", cell_id=" + str(self.cell_id) + ")"

        return repr_str

    def update(self, old_state, new_state):
        """Update state dictionary"""

        self.states[old_state] = self.states[old_state] - 1
        self.states[new_state] = self.states[new_state] + 1

def read_host_files(host_pos_files, init_cond_files, region_files, states, sim_type="INDIVIDUAL"):
    """Read all files associated with host initial state: position, initial state, region files."""

    if sim_type == "INDIVIDUAL":

        assert len(host_pos_files) == len(init_cond_files), "Number of input files do not match!"

        if region_files is not None:
            region_files = region_files.split(",")
            assert len(host_pos_files) == len(region_files), "Number of input files do not match!"
        else:
            region_files = [None for _ in host_pos_files]

        all_hosts = []

        for i, host_file in enumerate(host_pos_files):
            hosts = read_host_file(host_file, default_region=i, hostID_start=len(all_hosts))
            hosts = read_init_cond(hosts, init_cond_files[i])
            if region_files[i] is not None:
                read_regions(hosts, region_files[i])
            all_hosts.extend(hosts)

        return (all_hosts, None, None)

    if sim_type == "RASTER":

        # TODO implement raster region files

        assert len(host_pos_files) == 1, "Multiple host raster files provided!"

        assert len(init_cond_files) == 1, "Only provide file stub for initial condition rasters!"

        host_raster = raster_tools.RasterData.from_file(host_pos_files[0])
        state_rasters = []
        for state in states:
            state_rasters.append(raster_tools.RasterData.from_file(
                init_cond_files[0] + "_" + state + ".txt"))

        # Save raster header from host file
        header = host_raster.header_vals

        all_cells = []
        all_hosts = []
        cell_id = 0
        host_id = 0
        cellsize = host_raster.header_vals['cellsize']

        for row in range(host_raster.header_vals['nrows']):
            for col in range(host_raster.header_vals['ncols']):
                nhosts = host_raster.array[row, col]
                if nhosts <= 0:
                    continue
                hosts = []
                x = host_raster.header_vals['xllcorner'] + (col * cellsize) + (cellsize / 2)
                y = host_raster.header_vals['yllcorner'] + (
                    host_raster.header_vals['nrows'] * cellsize) - row * cellsize - (cellsize / 2)

                for i, state in enumerate(states):
                    nhosts_state = int(state_rasters[i].array[row, col])
                    for j in range(nhosts_state):
                        host = Host(x, y, state, host_id=host_id, cell_id=cell_id)
                        host_id += 1
                        hosts.append(host)
                        all_hosts.append(host)

                if len(hosts) != nhosts:
                    print(nhosts, len(hosts))
                    raise ValueError("Incorrect number of hosts!")

                all_cells.append(Cell((row, col), hosts, cell_id))
                cell_id += 1

        return (all_hosts, all_cells, header)

def read_sus_inf_files(all_cells, header, sus_file, inf_file, sim_type="INDIVIDUAL"):
    """Read all files associated with host susceptibility and infectiousness."""

    if sim_type == "INDIVIDUAL":

        return

    if sim_type == "RASTER":

        try:
            sus_raster = raster_tools.RasterData.from_file(sus_file)
        except FileNotFoundError:
            sus_raster = raster_tools.RasterData(
                shape=(header['nrows'], header['ncols']),
                llcorner=(header['xllcorner'], header['yllcorner']),
                cellsize=header['cellsize'],
                NODATA_value=header['NODATA_value'],
                array=np.ones((header['nrows'], header['ncols']))
            )

        try:
            inf_raster = raster_tools.RasterData.from_file(inf_file)
        except FileNotFoundError:
            inf_raster = raster_tools.RasterData(
                shape=(header['nrows'], header['ncols']),
                llcorner=(header['xllcorner'], header['yllcorner']),
                cellsize=header['cellsize'],
                NODATA_value=header['NODATA_value'],
                array=np.ones((header['nrows'], header['ncols']))
            )

        for cell in all_cells:
            row, col = cell.cell_position
            cell.susceptibility = sus_raster.array[row, col]
            cell.infectiousness = inf_raster.array[row, col]

def read_host_file(filename, default_region=0, hostID_start=0):
    """Read host file detailing host positions."""

    with open(filename, "r") as f:
        nhosts = int(f.readline())

        all_hosts = []

        for i in range(hostID_start, hostID_start + nhosts):
            x, y = f.readline().split()
            all_hosts.append(Host(float(x), float(y), host_id=i, reg=default_region))

    return all_hosts


def read_init_cond(all_hosts, filename):
    """Read host initial state file."""

    with open(filename, "r") as f:
        nhosts = int(f.readline())

        for i in range(nhosts):
            state = f.readline().strip()
            all_hosts[i].initialise_state(state)

    return all_hosts

def read_regions(all_hosts, filename):
    """Read region file giving regions to which hosts belong."""

    with open(filename, "r") as f:
        nhosts = int(f.readline())

        for i in range(nhosts):
            region = int(f.readline().strip())
            all_hosts[i].reg = region
