import importlib
import pandas as pd
import numpy as np


class InterventionHandler:

    def __init__(self, parent_sim):
        self.parent_sim = parent_sim
        self.interventions = []

        if self.parent_sim.params['InterventionScripts'] is not None:
            intervention_scripts = self.parent_sim.params['InterventionScripts'].split(",")
            intervention_freqs = self.parent_sim.params['InterventionUpdateFrequencies'].split(",")
            for script, freq in zip(intervention_scripts, intervention_freqs):
                intervention_module = importlib.import_module(script)
                intervention_module = importlib.reload(intervention_module)
                self.interventions.append(intervention_module.Intervention(
                    float(freq), self.parent_sim.params['init_hosts']))

            for i, intervention in enumerate(self.interventions):
                if intervention.type == "CONTINUOUS":
                    # Set up rate structures
                    self.parent_sim.rate_handler.add_rate_struct("Intervention_"+str(i),
                                                                 intervention.rate_size)

    @property
    def next_intervention_time(self):
        try:
            return min(self.next_interventions, default=np.inf)
        except AttributeError:
            raise AttributeError("Must first initialise_rates before accessing "
                                 "next_intervention_time!")

    def initialise_rates(self, all_hosts):
        """Set rates for continuous interventions, and setup next intervention times."""
        self.next_interventions = []

        for i, intervention in enumerate(self.interventions):
            if intervention.type == "CONTINUOUS":
                rates = intervention.update(all_hosts, 0)
                for j, rate in enumerate(rates):
                    self.parent_sim.rate_handler.insert_rate(j, rate, "Intervention_"+str(i))

            # Next intervention times
            self.next_interventions.append(intervention.update_freq)

    def update(self, all_hosts, time):
        """Carry out any intervention updates or actions for the specified time."""
        for i, int_time in enumerate(self.next_interventions):
            if time == int_time:
                if self.interventions[i].type == "CONTINUOUS":
                    rates = self.interventions[i].update(all_hosts, time)
                    for j, rate in enumerate(rates):
                        self.parent_sim.rate_handler.insert_rate(j, rate, "Intervention_"+str(i))
                else:
                    hostID, event_type = self.interventions[i].action(all_hosts, time)
                    self.parent_sim.event_handler.do_event(event_type, hostID, all_hosts)
                self.next_interventions[i] += self.interventions[i].update_freq

    def update_on_event(self, all_hosts, time):
        """Carry out updates after an event is carried out."""
        for i, intervention in enumerate(self.interventions):
            if intervention.type == "CONTINUOUS":
                rates = intervention.update(all_hosts, time)
                for j, rate in enumerate(rates):
                    self.parent_sim.rate_handler.insert_rate(j, rate, "Intervention_"+str(i))
            else:
                intervention.update(all_hosts, time)

    def action(self, event_type, eventID, all_hosts, all_cells):
        """Carry out an action for a continuous intervention."""

        int_num = int(event_type.split("_")[1])
        time = self.parent_sim.time

        hostID, event_type = self.interventions[int_num].action(all_hosts, time, eventID)
        self.parent_sim.event_handler.do_event(event_type, hostID, all_hosts)

    def output(self, iteration=0, file_stub="output"):
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
