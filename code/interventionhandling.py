"""Module for handling control interventions in Individual Epidemic Simulator."""

import pdb
import importlib
import pandas as pd
import numpy as np


class InterventionHandler:
    """Intervention handling class to keep interventions up to date and carry out events."""

    def __init__(self, parent_sim):
        self.parent_sim = parent_sim
        self.interventions = []
        self.next_interventions = None

        scripts = self.parent_sim.params['InterventionScripts']

        if scripts is not None:
            intervention_freqs = self.parent_sim.params['InterventionUpdateFrequencies']
            intervention_options = self.parent_sim.params['InterventionOptions']
            if isinstance(intervention_options, str):
                intervention_options = intervention_options.split(',')
            elif isinstance(intervention_options, list):
                pass
            else:
                raise TypeError("InterventionOptions must be list or string!")

            if isinstance(scripts, str):
                intervention_scripts = scripts.split(",")
                if intervention_freqs is None:
                    intervention_freqs = [None]*len(intervention_scripts)
                else:
                    intervention_freqs = [float(x) for x in intervention_freqs.split(",")]
                zipped_data = zip(intervention_scripts, intervention_freqs, intervention_options)
                for script, freq, options in zipped_data:
                    intervention_module = importlib.import_module(script)
                    intervention_module = importlib.reload(intervention_module)
                    self.interventions.append(intervention_module.Intervention(
                        freq, self.parent_sim.params['init_hosts'],
                        self.parent_sim.params['init_cells']), options)

            elif isinstance(scripts, list):
                if intervention_freqs is None:
                    intervention_freqs = [None]*len(scripts)
                else:
                    intervention_freqs = [float(x) for x in intervention_freqs.split(",")]
                zipped_data = zip(scripts, intervention_freqs, intervention_options)
                for intervention_class, freq, options in zipped_data:
                    self.interventions.append(
                        intervention_class(
                            freq, self.parent_sim.params['init_hosts'],
                            self.parent_sim.params['init_cells'], options))

            else:
                raise ValueError("InterventionScripts must be list or string!")

            for i, intervention in enumerate(self.interventions):
                if intervention.type == "CONTINUOUS":
                    # Set up rate structures
                    rate_factor = getattr(intervention, "rate_factor", None)
                    self.parent_sim.rate_handler.add_rate_struct("Intervention_"+str(i),
                                                                 intervention.rate_size,
                                                                 rate_factor=rate_factor)

    @property
    def next_intervention_time(self):
        """Time of next intervention update for each intervention."""

        try:
            return min([x for x in self.next_interventions if x is not None], default=np.inf)
        except TypeError:
            raise TypeError("Must first initialise_rates before accessing "
                            "next_intervention_time!")

    def initialise_rates(self, all_hosts, all_cells):
        """Set rates for continuous interventions, and setup next intervention times."""
        self.next_interventions = []

        for i, intervention in enumerate(self.interventions):
            if intervention.type == "CONTINUOUS":
                rate_updates = intervention.update(all_hosts, 0, all_cells, initial=True)
                for j, rate in rate_updates:
                    self.parent_sim.rate_handler.insert_rate(j, rate, "Intervention_"+str(i))

            # Next intervention times
            self.next_interventions.append(intervention.update_freq)

    def update(self, all_hosts, time, all_cells):
        """Carry out any intervention updates or actions for the specified time.

        This function will be called every update frequency.
        """

        for i, int_time in enumerate(self.next_interventions):
            if time == int_time:
                if self.interventions[i].type == "CONTINUOUS":
                    get_rate_fn = lambda x: self.parent_sim.rate_handler.get_rate(
                        x, "Intervention_"+str(i))
                    rate_updates = self.interventions[i].update(all_hosts, time, all_cells,
                                                                get_rate_fn=get_rate_fn)
                    for rate_id, rate in rate_updates:
                        self.parent_sim.rate_handler.insert_rate(
                            rate_id, rate, "Intervention_"+str(i))
                else:
                    events = self.interventions[i].update(all_hosts, time, all_cells)
                    for host_id, event_type in events:
                        self.parent_sim.event_handler.do_event(
                            event_type, host_id, all_hosts, all_cells)
                self.next_interventions[i] += self.interventions[i].update_freq

    def update_on_event(self, event, all_hosts, time, all_cells):
        """Carry out updates after an event is carried out."""
        for i, intervention in enumerate(self.interventions):
            if intervention.type == "CONTINUOUS":
                rate_fn = lambda x: self.parent_sim.rate_handler.get_rate(
                    x, "Intervention_"+str(i))
                rate_updates = intervention.update(all_hosts, time, all_cells, after_event=event,
                                                   get_rate_fn=rate_fn)
                for rate_id, rate in rate_updates:
                    self.parent_sim.rate_handler.insert_rate(rate_id, rate, "Intervention_"+str(i))
            else:
                events = intervention.update(all_hosts, time, all_cells, after_event=event)
                for host_id, event_type in events:
                    self.parent_sim.event_handler.do_event(
                        event_type, host_id, all_hosts, all_cells)

    def action(self, intervention_id, event_id, all_hosts, all_cells):
        """Carry out an action for a continuous intervention."""

        int_num = int(intervention_id.split("_")[1])
        time = self.parent_sim.time

        events = self.interventions[int_num].action(all_hosts, time, event_id, all_cells)
        for host_id, event_type in events[:-1]:
            event_done = self.parent_sim.event_handler.do_event(
                event_type, host_id, all_hosts, all_cells)
            self.update_on_event(event_done, all_hosts, time, all_cells)

        event_done = self.parent_sim.event_handler.do_event(
            events[-1][1], events[-1][0], all_hosts, all_cells)

        return event_done

    def output(self, iteration=0, file_stub="output"):
        """Save any output generated by interventions."""

        intervention_data = []

        for i, intervention in enumerate(self.interventions):
            try:
                intervention_data.append(intervention.output())
            except AttributeError:
                intervention_data.append(None)

            filename = file_stub + "_Intervention" + str(i) + "_" + str(iteration) + ".csv"
            if isinstance(intervention_data[-1], pd.DataFrame):
                if self.parent_sim.params['OutputFiles'] is True:
                    intervention_data[-1].to_csv(filename, index=False)
            elif intervention_data[-1] is not None:
                raise TypeError("Intervention output format must be pandas DataFrame!")

        if len(intervention_data) == 1:
            intervention_data = intervention_data[0]

        return intervention_data
