"""Base intervention class for continuous removal from regions, subject to a budget constraint.

Behaviour can be adjusted by changing the config_params using the intervention_utilties within the
utilities module.
"""

import numpy as np
import pandas as pd

config_params = {'budget': 0, 'nregions': 3, 'priorities': [2, 1, 0], 'removal_rate': 1.0}


class Intervention:
    def __init__(self, update_freq, all_hosts):
        # Set update frequency
        self.update_freq = update_freq

        # Set intervention type
        self.type = "CONTINUOUS"

        # Set size of required rate structure
        self.rate_size = config_params['nregions']

        # Create region -> to hostID map
        self.region_map = {key: [] for key in range(config_params['nregions'])}
        for host in all_hosts:
            region = host.reg
            self.region_map[region].append(host.hostID)

        # Create log of expenditure
        self.log = []

    def action(self, all_hosts, time, eventID):
        """Carry out event."""
        # select random infected host in correct region, and carry out cull event
        nhosts_region = len(self.region_map[eventID])
        indx = np.random.randint(0, nhosts_region)
        host = all_hosts[self.region_map[eventID][indx]]
        while host.state != "I":
            indx = np.random.randint(0, nhosts_region)
            host = all_hosts[self.region_map[eventID][indx]]

        return (host.hostID, "Cull")

    def update(self, all_hosts, time):
        """Adjust rates."""
        # Adjust rates to account for changing budget and/or changing number of infected hosts

        already_spent = 0.0
        rates = [0]*config_params['nregions']
        expenses = [0]*config_params['nregions']

        for region in config_params['priorities']:
            N_inf = len([i for i in self.region_map[region] if all_hosts[i].state == "I"])
            expenses[region] = np.max([0, np.min([
                N_inf, config_params['budget'] - already_spent])])
            already_spent += expenses[region]
            rates[region] = expenses[region]*config_params['removal_rate']

        if time == 0:
            # Start of new epidemic, reset log
            self.log = [[time] + expenses]
        elif expenses != self.log[-1][1:]:
            # If expenses have changed, save to log
            self.log.append([time] + expenses)

        return rates

    def output(self):
        """Output control log to data frame."""

        data = pd.DataFrame(
            self.log,
            columns=['time'] + ['ExpenseRegion'+str(region) for region in
                                range(config_params['nregions'])])

        return data
