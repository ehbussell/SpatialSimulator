"""Tools for visualising input and output from the simulator."""

import time
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import colors
from matplotlib.animation import FuncAnimation
from ..code import hosts

plt.style.use("seaborn-whitegrid")
graphParams = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 8,
   'xtick.labelsize': 8,
   'ytick.labelsize': 8,
   'text.usetex': False,
   'figure.figsize': [6, 4]
}


def plot_host_landscape(host_file, region_file=None, save_file=None):
    """Generate plot of host position file.

    Arguments:
        host_file:      Host position filename.  Alternatively, can specify a list of filenames for
                        hosts contained withing multiple files.  Note: if multiple filenames are
                        specified then each file is assumed to be a separate region, unless a
                        corresponding list of region files is given.
        region_file:    Filename for region information on hosts.  Can specify a list of region
                        files if a corresponding list of host files is given.
        save_file:      If specified, plot is saved to this filename.
    """

    all_hosts = hosts.read_host_file(host_file)

    if region_file is not None:
        hosts.read_regions(all_hosts, region_file)

    regions = set(host.reg for host in all_hosts)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for region in regions:
        x, y = zip(*[(host.xpos, host.ypos) for host in all_hosts if host.reg == region])
        ax1.plot(x, y, '.', label="Region " + str(region))
    ax1.set_xlabel("X Position")
    ax1.set_ylabel("Y Position")

    ax1.set_aspect("equal")
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax1.legend(loc="center left", bbox_to_anchor=(1, 0.5), fancybox=True)

    ymin, ymax = ax1.get_ylim()
    ax1.set_ylim([ymin - 0.05*(ymax-ymin), ymax + 0.05*(ymax-ymin)])
    xmin, xmax = ax1.get_xlim()
    ax1.set_xlim([xmin - 0.05*(xmax-xmin), xmax + 0.05*(xmax-xmin)])

    fig.tight_layout()

    if save_file is not None:
        fig.savefig(save_file, dpi=300)

    plt.show()


def plot_results(hosts_filename="output_hosts_0.csv", event_filename="output_events_0.csv"):
    """Animate a particular epidemic simulation."""

    state_colours = {
        "S": colors.to_rgba("green"),
        "E": colors.to_rgba("cyan"),
        "C": colors.to_rgba("blue"),
        "D": colors.to_rgba("orange"),
        "I": colors.to_rgba("red"),
        "R": colors.to_rgba("grey"),
        "Culled": colors.to_rgba("black"),
    }

    hosts_df = pd.read_csv(hosts_filename)
    events_df = pd.read_csv(event_filename)

    for i, x in enumerate(hosts_df.columns):
        if "timeEnter" in x:
            host_start_idx = i
            break

    host_points = np.zeros(len(hosts_df), dtype=[('position', float, 2),
                                                 ('colour', float, 4)])

    for index, host in hosts_df.iterrows():
        host_points['position'][index] = (host['posX'], host['posY'])
        for name, timeEnter in host[host_start_idx:].iteritems():
            if 0 in eval(timeEnter):
                state = name[9:]
        host_points['colour'][index] = state_colours[state]

    initial_colours = host_points['colour'][:]

    frame_rate = 30
    animation_length = 5
    frame_interval = max(events_df['time']) / (animation_length * frame_rate)
    nframes = animation_length * frame_rate + 1

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax.set_xlim(min(host_points['position'][:, 0]), max(host_points['position'][:, 0]))
    ax.set_xticks([])
    ax.set_ylim(min(host_points['position'][:, 1]), max(host_points['position'][:, 1]))
    ax.set_yticks([])

    scat = ax.scatter(host_points['position'][:, 0], host_points['position'][:, 1],
                      c=host_points['colour'])

    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, weight="bold", fontsize=12,
                        bbox=dict(facecolor='white', alpha=0.6))

    global start_from_idx
    start_from_idx = 0

    def init():
        global start_from_idx
        start_from_idx = 0
        host_points['colour'] = initial_colours
        new_time = 0
        scat.set_color(host_points['colour'])
        time_text.set_text('time = %.3f' % new_time)

        return scat, time_text

    def update(frame_number):
        new_time = frame_number * frame_interval
        global start_from_idx

        row = events_df.iloc[start_from_idx]

        while row['time'] <= new_time:
            hostID = row['hostID']
            new_state = row['newState']
            host_points['colour'][hostID] = state_colours[new_state]
            start_from_idx += 1
            if start_from_idx >= len(events_df):
                break
            row = events_df.iloc[start_from_idx]

        scat.set_color(host_points['colour'])
        time_text.set_text('time = %.3f' % new_time)

        return scat, time_text

    animation = FuncAnimation(fig, update, interval=1000/frame_rate, frames=nframes, blit=True,
                              repeat=False, init_func=init)
    plt.show()
