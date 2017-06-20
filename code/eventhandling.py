import numpy as np


class EventHandler:

    def __init__(self, parent_sim, rate_handler):
        self.parent_sim = parent_sim
        self.rate_handler = rate_handler
        self.cache_kernel = self.parent_sim.params["CacheKernel"]

        if self.cache_kernel is True:
            self.kernel = self.kernel_cached
        else:
            self.kernel = self.kernel_uncached

    def kernel_cached(self, i, hostID):
        return self.parent_sim.params['kernel_vals'][i, hostID]

    def kernel_uncached(self, i, hostID):
        return self.parent_sim.params['kernel'](np.linalg.norm(
            [all_hosts[i].x-all_hosts[hostID].x,
             all_hosts[i].y-all_hosts[hostID].y]))

    def do_event(self, event_type, eventID, all_hosts):
        if event_type == "Infection" or event_type == "Advance":
            self.do_event_standard(eventID, all_hosts)
        elif event_type == "Cull":
            self.do_event_cull(eventID, all_hosts)
        elif event_type.startswith("Intervention"):
            self.parent_sim.intervention_handler.action(event_type, eventID, all_hosts)
        else:
            raise ValueError("Unrecognised event type: {}!".format(event_type))

    def do_event_standard(self, hostID, all_hosts):
        old_state = all_hosts[hostID].state
        new_state = self.parent_sim.params['next_state'](old_state)
        all_hosts[hostID].state = new_state

        if old_state == "S":
            self.rate_handler.insert_rate(hostID, 0.0, "Infection")
        if new_state in "ECDI":
            self.rate_handler.insert_rate(
                hostID, self.parent_sim.params[new_state + 'AdvRate'], "Advance")
        if new_state == "R":
            self.rate_handler.insert_rate(hostID, 0.0, "Advance")
            # Distribute rate changes
            for i in range(len(all_hosts)):
                if all_hosts[i].state == "S":
                    old_rate = self.rate_handler.get_rate(i, "Infection")
                    new_rate = old_rate - self.kernel(i, hostID)
                    self.rate_handler.insert_rate(i, new_rate, "Infection")

        if new_state in "CI":
            # Distribute rate changes
            for i in range(len(all_hosts)):
                if all_hosts[i].state == "S":
                    old_rate = self.rate_handler.get_rate(i, "Infection")
                    new_rate = old_rate + self.kernel(i, hostID)
                    self.rate_handler.insert_rate(i, new_rate, "Infection")

        all_hosts[hostID].jump_times[new_state] = self.parent_sim.time

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, hostID, old_state, new_state))
        self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][old_state] -= 1
        self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][new_state] += 1

    def do_event_cull(self, hostID, all_hosts):
        old_state = all_hosts[hostID].state
        new_state = "Culled"
        all_hosts[hostID].state = new_state

        self.rate_handler.insert_rate(hostID, 0.0, "Infection")
        self.rate_handler.insert_rate(hostID, 0.0, "Advance")
        if old_state in "CI":
            # Distribute rate changes
            for i in range(len(all_hosts)):
                if all_hosts[i].state == "S":
                    old_rate = self.rate_handler.get_rate(i, "Infection")
                    new_rate = old_rate - self.kernel(i, hostID)
                    self.rate_handler.insert_rate(i, new_rate, "Infection")

        all_hosts[hostID].jump_times[new_state] = self.parent_sim.time

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, hostID, old_state, new_state))
        self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][old_state] -= 1
        self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][new_state] += 1
