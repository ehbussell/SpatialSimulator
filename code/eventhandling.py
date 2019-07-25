"""Module for handling events in Individual Epidemic Simulator."""

import pdb
import numpy as np

class EventHandler:
    """Class to carry out all events on hosts and cells."""

    def __init__(self, parent_sim, rate_handler):
        self.parent_sim = parent_sim
        self.rate_handler = rate_handler
        self.cache_kernel = self.parent_sim.params["CacheKernel"]

        if self.parent_sim.params['SimulationType'] == "INDIVIDUAL":
            self.do_event_advance = self.do_event_standard
            self.do_event_infection = self.do_event_standard

            if self.cache_kernel is True:
                self.kernel = self.kernel_cached
            else:
                self.kernel = self.kernel_uncached

        elif self.parent_sim.params['SimulationType'] == "RASTER":
            self.do_event_advance = self.do_event_adv_raster
            self.do_event_infection = self.do_event_inf_raster

            self.kernel = self.kernel_raster

        else:
            raise ValueError("Unrecognised SimulationType!")

    def kernel_cached(self, host_id1, host_id2):
        """Calculate kernel between two hosts - from cached value."""

        return self.parent_sim.params['kernel_vals'][host_id1, host_id2]

    def kernel_uncached(self, host_id1, host_id2):
        """Calculate kernel between two hosts - calculated on-they-fly."""

        # TODO FIXME uncached kernels
        pass
        # return self.parent_sim.params['kernel'](np.linalg.norm(
        #     [all_hosts[i].x-all_hosts[hostID].x,
        #      all_hosts[i].y-all_hosts[hostID].y]))

    def kernel_raster(self, cell_rel_pos):
        """Get kernel value in raster mode for given relative position."""

        shape = self.parent_sim.params['coupled_kernel'].shape

        return self.parent_sim.params['coupled_kernel'][cell_rel_pos[0] + int(shape[0]/2),
                                                        cell_rel_pos[1] + int(shape[1]/2)]

    def do_event(self, event_type, event_id, all_hosts, all_cells):
        """Carry out event on given host or cell."""

        if event_type == "Infection":
            return self.do_event_infection(event_id, all_hosts, all_cells)
        elif event_type == "Advance":
            return self.do_event_advance(event_id, all_hosts, all_cells)
        elif event_type == "Sporulation":
            return self.do_event_sporulation(event_id, all_hosts, all_cells)
        elif event_type == "CULL":
            return self.do_event_cull(event_id, all_hosts, all_cells)
        elif event_type.startswith("Intervention"):
            return self.parent_sim.intervention_handler.action(event_type, event_id, all_hosts,
                                                               all_cells)
        else:
            raise ValueError("Unrecognised event type: {}!".format(event_type))

    def do_event_standard(self, host_id, all_hosts, all_cells):
        """Carry out standard infection/advance event in Individual model mode."""

        old_state = all_hosts[host_id].state
        new_state = self.parent_sim.params['next_state'](old_state)

        if old_state == "S":
            self.rate_handler.insert_rate(host_id, 0.0, "Infection")
        if new_state in "ECDI":
            self.rate_handler.insert_rate(
                host_id, self.parent_sim.params[new_state + 'AdvRate'], "Advance")
        if new_state == "R":
            self.rate_handler.insert_rate(host_id, 0.0, "Advance")
            # Distribute rate changes
            self.distribute_removal_individual(host_id, all_hosts)
        if new_state in "CI":
            # Distribute rate changes
            self.distribute_infection_individual(host_id, all_hosts)

        all_hosts[host_id].update_state(new_state, self.parent_sim.time)

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, host_id, old_state, new_state))
        # self.parent_sim.run_params['region_summary'][all_hosts[host_id].reg][old_state] -= 1
        # self.parent_sim.run_params['region_summary'][all_hosts[host_id].reg][new_state] += 1

        return (host_id, None, old_state, new_state)

    def do_event_adv_raster(self, host_id, all_hosts, all_cells):
        """Carry out advance event in raster model."""

        old_state = all_hosts[host_id].state
        new_state = self.parent_sim.params['next_state'](old_state)

        cell = all_cells[all_hosts[host_id].cell_id]

        all_hosts[host_id].update_state(new_state, self.parent_sim.time)
        cell.update(old_state, new_state)

        if new_state in "ECDI":
            self.rate_handler.insert_rate(
                host_id, self.parent_sim.params[new_state + 'AdvRate'], "Advance")
        if new_state == "R":
            self.rate_handler.insert_rate(host_id, 0.0, "Advance")
            # Distribute rate changes to coupled cells
            self.distribute_removal_raster(host_id, all_hosts, all_cells)

        if new_state in "CI":
            # Distribute rate changes to coupled cells
            self.distribute_infection_raster(host_id, all_hosts, all_cells)

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, host_id, old_state, new_state))

        return (host_id, cell, old_state, new_state)

    def do_event_inf_raster(self, cell_id, all_hosts, all_cells):
        """Carry out infection event in Raster model."""

        cell = all_cells[cell_id]

        nsus = cell.states["S"]
        if nsus <= 0:
            raise ValueError("No susceptibles to infect!\n" + str(cell.states))

        for host in cell.hosts:
            if host.state == "S":
                host_id = host.host_id
                break

        new_state = self.parent_sim.params['next_state']("S")
        all_hosts[host_id].update_state(new_state, self.parent_sim.time)
        cell.update("S", new_state)

        # TODO handle if sporulation starts within same cell
        # TODO combine this into distribution of rates
        old_inf_rate = self.rate_handler.get_rate(cell_id, "Infection")
        new_inf_rate = old_inf_rate * ((nsus - 1) / nsus)
        self.rate_handler.insert_rate(cell_id, new_inf_rate, "Infection")

        if new_state in "ECDI":
            self.rate_handler.insert_rate(
                host_id, self.parent_sim.params[new_state + 'AdvRate'], "Advance")

        if new_state in "CI":
            # Distribute rate changes to coupled cells
            self.distribute_infection_raster(host_id, all_hosts, all_cells)

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, host_id, "S", new_state))

        return (host_id, cell_id, "S", new_state)

    def do_event_sporulation(self, cell_id, all_hosts, all_cells, debug=False):
        """Carry out sporulation event in raster model."""

        selection_val = np.random.random_sample()
        selected_idx = self.parent_sim.params['vs_kernel'].select_event(selection_val)

        kernel_shape = self.parent_sim.params['kernel'].shape
        kernel_pos = np.unravel_index(selected_idx, kernel_shape)
        cell_rel_pos = [kernel_pos[i] - int(kernel_shape[i]/2) for i in range(2)]

        random_num = np.random.random_sample()

        cell_pos = tuple(item1 + item2 for item1, item2
                         in zip(all_cells[cell_id].cell_position, cell_rel_pos))

        cell_id = self.parent_sim.params['cell_map'].get(cell_pos, None)

        if cell_id is not None:
            random_num = np.random.random_sample()
            if random_num < (all_cells[cell_id].states["S"] / self.parent_sim.params['MaxHosts']):
                event = self.do_event("Infection", cell_id, all_hosts, all_cells)
                return event

        if debug:
            return cell_id

    def do_event_cull(self, host_id, all_hosts, all_cells):
        """Carry out cull intervention, removing a host."""

        old_state = all_hosts[host_id].state
        new_state = "Culled"
        cell_id = all_hosts[host_id].cell_id
        all_hosts[host_id].update_state(new_state, self.parent_sim.time)

        self.rate_handler.insert_rate(host_id, 0.0, "Advance")

        if cell_id is None:
            self.rate_handler.insert_rate(host_id, 0.0, "Infection")
        else:
            all_cells[cell_id].update(old_state, new_state)

        if old_state in "CI":
            # Distribute rate changes
            if cell_id is None:
                self.distribute_removal_individual(host_id, all_hosts)
            else:
                self.distribute_removal_raster(host_id, all_hosts, all_cells)

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, host_id, old_state, new_state))
        # self.parent_sim.run_params['region_summary'][all_hosts[host_id].reg][old_state] -= 1
        # self.parent_sim.run_params['region_summary'][all_hosts[host_id].reg][new_state] += 1


        return (host_id, cell_id, old_state, new_state)

    def distribute_infection_individual(self, host_id, all_hosts):
        """Host has just become infectious - distribute rate changes in Individual model."""

        for host in all_hosts:
            if host.state == "S":
                coupled_host_id = host.host_id
                old_rate = self.rate_handler.get_rate(coupled_host_id, "Infection")
                new_rate = old_rate + self.kernel(coupled_host_id, host_id)
                self.rate_handler.insert_rate(coupled_host_id, new_rate, "Infection")

    def distribute_infection_raster(self, host_id, all_hosts, all_cells):
        """Host has just become infectious - distribute rate changes in Raster model."""

        cell = all_cells[all_hosts[host_id].cell_id]

        # Update sporulation rate
        if self.parent_sim.params['VirtualSporulationStart'] is not None:
            self.rate_handler.insert_rate(
                cell.cell_id,
                cell.states["C"] + cell.states["I"], "Sporulation")

        # Update coupled cells
        for cell2_rel_pos in self.parent_sim.params['coupled_positions']:
            cell2_pos = tuple(item1 + item2 for item1, item2
                              in zip(cell.cell_position, cell2_rel_pos))
            cell2_id = self.parent_sim.params['cell_map'].get(cell2_pos, None)
            if cell2_id is None:
                continue
            cell2 = all_cells[cell2_id]
            if cell2.states["S"] > 0:
                old_rate = self.rate_handler.get_rate(cell2_id, "Infection")
                new_rate = old_rate + (self.kernel(cell2_rel_pos) *
                                       cell2.states["S"]) / self.parent_sim.params['MaxHosts']
                self.rate_handler.insert_rate(cell2_id, new_rate, "Infection")

    def distribute_removal_individual(self, host_id, all_hosts):
        """Host has just lost infectivity - distribute rate changes in Individual model."""

        for host in all_hosts:
            if host.state == "S":
                coupled_host_id = host.host_id
                old_rate = self.rate_handler.get_rate(coupled_host_id, "Infection")
                new_rate = old_rate - self.kernel(coupled_host_id, host_id)
                self.rate_handler.insert_rate(coupled_host_id, new_rate, "Infection")

    def distribute_removal_raster(self, host_id, all_hosts, all_cells):
        """Host has just lost infectivity - distribute rate changes in Raster model."""

        cell = all_cells[all_hosts[host_id].cell_id]

        # Update sporulation rate
        if self.parent_sim.params['VirtualSporulationStart'] is not None:
            self.rate_handler.insert_rate(
                cell.cell_id,
                cell.states["C"] + cell.states["I"], "Sporulation")

        # Update coupled cells
        for cell2_rel_pos in self.parent_sim.params['coupled_positions']:
            cell2_pos = tuple(item1 + item2 for item1, item2
                              in zip(cell.cell_position, cell2_rel_pos))
            cell2_id = self.parent_sim.params['cell_map'].get(cell2_pos, None)
            if cell2_id is None:
                continue
            cell2 = all_cells[cell2_id]
            if cell2.states["S"] > 0:
                old_rate = self.rate_handler.get_rate(cell2_id, "Infection")
                new_rate = old_rate - (self.kernel(cell2_rel_pos) *
                                       cell2.states["S"])
                self.rate_handler.insert_rate(cell2_id, new_rate, "Infection")
