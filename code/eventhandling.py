import numpy as np


class EventHandler:

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

    def kernel_cached(self, i, hostID):
        return self.parent_sim.params['kernel_vals'][i, hostID]

    def kernel_uncached(self, i, hostID):
        # TODO fix uncached kernels
        pass
        # return self.parent_sim.params['kernel'](np.linalg.norm(
        #     [all_hosts[i].x-all_hosts[hostID].x,
        #      all_hosts[i].y-all_hosts[hostID].y]))

    def kernel_raster(self, cell1_pos, cell2_pos):
        # row_diff = abs(cell1_pos[0] - cell2_pos[0])
        # col_diff = abs(cell1_pos[1] - cell2_pos[1])

        return self.parent_sim.params['kernel'].array[abs(cell1_pos[0] - cell2_pos[0]),
                                                      abs(cell1_pos[1] - cell2_pos[1])]

    def do_event(self, event_type, eventID, all_hosts, all_cells):
        if event_type == "Infection":
            self.do_event_infection(eventID, all_hosts, all_cells)
        elif event_type == "Advance":
            self.do_event_advance(eventID, all_hosts, all_cells)
        elif event_type == "Cull":
            self.do_event_cull(eventID, all_hosts, all_cells)
        elif event_type.startswith("Intervention"):
            self.parent_sim.intervention_handler.action(event_type, eventID, all_hosts, all_cells)
        else:
            raise ValueError("Unrecognised event type: {}!".format(event_type))

    def do_event_standard(self, hostID, all_hosts, all_cells):
        old_state = all_hosts[hostID].state
        new_state = self.parent_sim.params['next_state'](old_state)

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

        all_hosts[hostID].update_state(new_state, self.parent_sim.time)

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, hostID, old_state, new_state))
        # self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][old_state] -= 1
        # self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][new_state] += 1

    def do_event_adv_raster(self, hostID, all_hosts, all_cells):
        old_state = all_hosts[hostID].state
        new_state = self.parent_sim.params['next_state'](old_state)

        all_hosts[hostID].update_state(new_state, self.parent_sim.time)
        all_cells[all_hosts[hostID].cell_id].update(old_state, new_state)

        if new_state in "ECDI":
            self.rate_handler.insert_rate(
                hostID, self.parent_sim.params[new_state + 'AdvRate'], "Advance")
        if new_state == "R":
            self.rate_handler.insert_rate(hostID, 0.0, "Advance")
            # Distribute rate changes
            for cell_id, cell in enumerate(all_cells):
                if cell.states.get("S", 0) > 0:
                    old_rate = self.rate_handler.get_rate(cell_id, "Infection")
                    new_rate = old_rate - (self.kernel(cell.cell_position,
                                           all_cells[all_hosts[hostID].cell_id].cell_position) *
                                           cell.states.get("S", 0))
                    self.rate_handler.insert_rate(cell_id, new_rate, "Infection")

        if new_state in "CI":
            # Distribute rate changes
            for cell_id, cell in enumerate(all_cells):
                if cell.states.get("S", 0) > 0:
                    old_rate = self.rate_handler.get_rate(cell_id, "Infection")
                    new_rate = old_rate + (self.kernel(
                        cell.cell_position, all_cells[all_hosts[hostID].cell_id].cell_position) * 
                                           cell.states.get("S", 0))
                    self.rate_handler.insert_rate(cell_id, new_rate, "Infection")

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, hostID, old_state, new_state))

    def do_event_inf_raster(self, cellID, all_hosts, all_cells):
        nsus = all_cells[cellID].states.get("S", 0)
        if nsus <= 0:
            raise ValueError("No susceptibles to infect!")

        for host in all_cells[cellID].hosts:
            if host.state == "S":
                hostID = host.host_id
                break
        
        new_state = self.parent_sim.params['next_state']("S")
        all_hosts[hostID].update_state(new_state, self.parent_sim.time)
        all_cells[cellID].update("S", new_state)

        old_inf_rate = self.rate_handler.get_rate(cellID, "Infection")
        new_inf_rate = old_inf_rate * ((nsus - 1) / nsus)
        self.rate_handler.insert_rate(cellID, new_inf_rate, "Infection")

        if new_state in "ECDI":
            self.rate_handler.insert_rate(
                hostID, self.parent_sim.params[new_state + 'AdvRate'], "Advance")

        if new_state in "CI":
            # Distribute rate changes
            for cell_id, cell in enumerate(all_cells):
                if cell.states.get("S", 0) > 0:
                    old_rate = self.rate_handler.get_rate(cell_id, "Infection")
                    new_rate = old_rate + (
                        self.kernel(cell.cell_position, all_cells[cellID].cell_position) *
                        cell.states.get("S", 0))
                    self.rate_handler.insert_rate(cell_id, new_rate, "Infection")

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, hostID, "S", new_state))

    def do_event_cull(self, hostID, all_hosts, all_cells):
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

        all_hosts[hostID].update_state(new_state, self.parent_sim.time)

        self.parent_sim.run_params['all_events'].append(
            (self.parent_sim.time, hostID, old_state, new_state))
        # self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][old_state] -= 1
        # self.parent_sim.run_params['region_summary'][all_hosts[hostID].reg][new_state] += 1
