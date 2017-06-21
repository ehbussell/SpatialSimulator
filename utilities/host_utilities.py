"""Tools for creating and editing host landscape files.

Methods:
    - create_rand_landscape
    - create_rand_initial_conditions
    - create_regions
"""

import itertools
import numpy as np


def create_rand_landscape(filename, nhosts, x_limits=[0, 1], y_limits=[0, 1]):
    """Create a host landscape with uniformly distributed host positions.

    Arguments:
        filename:   Name of file to which host data will be printed
        nhosts:     The number of hosts to randomly position
        x_limits:   Hosts are positioned between these limits (length 2 list) in x direction
        x_limits:   Hosts are positioned between these limits (length 2 list) in x direction
    """

    all_x = x_limits[0] + np.random.random_sample(nhosts)*x_limits[1]
    all_y = y_limits[0] + np.random.random_sample(nhosts)*y_limits[1]

    with open(filename, "w") as f:
        f.write(str(nhosts) + "\n")
        for i in range(nhosts):
            f.write(str(all_x[i]) + " " + str(all_y[i]) + "\n")


def create_grid_landscape(filename, host_layout=[10, 10], x_spacing=1, y_spacing=1,
                          llcorner=[0, 0]):
    """Create a host landscape with hosts positioned on a grid.

    Arguments:
        filename:       Name of file to which host data will be printed
        host_layout:    List giving dimensions of host grid [x, y].  Rows will have x hosts on
                        and columns y hosts
        x_spacing:      Distance between hosts in x direction
        y_spacing:      Distance between hosts in y direction
        llcorner:       [x, y] postion of lower left corner of grid
    """

    x_vals = llcorner[0] + np.array(range(host_layout[0])) * x_spacing
    y_vals = llcorner[1] + np.array(range(host_layout[1])) * y_spacing
    nhosts = host_layout[0] * host_layout[1]

    with open(filename, "w") as f:
        f.write(str(nhosts) + "\n")
        for y_val in y_vals:
            for x_val in x_vals:
                f.write(str(x_val) + " " + str(y_val) + "\n")


def create_rand_initial_conditions(filename, nhosts, rand_infs=0, fixed_infs=None):
    """Create a initial conditions files, detailing the initial state for each host.

    Hosts are either susceptible (S) or infected (I).  Hosts can be made infectious by random
    choice, fixed indices, or a combination of the two.

    Arguments:
        filename:   Name of file to which host data will be printed
        nhosts:     The number of hosts to assign an initial state
        rand_infs:  The number of hosts to randomly select to start infected
        fixed_infs: A list of indices for the hosts which will start infected
    """

    host_list = list(range(nhosts))
    infected = []

    if fixed_infs is not None:
        infected = fixed_infs
        host_list = [x for x in host_list if x not in infected]

    infected += list(np.random.choice(host_list, rand_infs, replace=False))

    with open(filename, "w") as f:
        f.write(str(nhosts) + "\n")
        for i in range(nhosts):
            if i in infected:
                f.write("I\n")
            else:
                f.write("S\n")


def create_regions(filename, nhosts, nregions, host_regions):
    """Create a region description file, detailing the containing region for each host.

    Arguments:
        filename:       Name of file to which host data will be printed
        nhosts:         The number of hosts to assign to a region
        nregions:       The total number of regions
        host_regions:   Either a list of regions for each host, or a list of number of hosts in
                        each region.  If the latter, then hosts are assigned to the regions in
                        index order.

    """

    if len(host_regions) == nregions:
        host_regions = list(itertools.chain.from_iterable(
            [[i]*val for i, val in enumerate(host_regions)]))
    elif len(host_regions) != nhosts:
        raise ValueError("Unknown host_regions length!")

    with open(filename, "w") as f:
        f.write(str(nhosts) + "\n")
        for i in range(nhosts):
                f.write(str(host_regions[i]) + "\n")
