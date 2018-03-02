"""Template for Intervention Class script."""

# Insert any module level code

# Must define an Intervention class which is initialised and used by the simulator
class Intervention:
    """Class to hold intervention methods."""

    # Initialisation sets up defining features of intervention
    def __init__(self, update_freq, all_hosts, all_cells=None):
        # Initialisation takes the update frequency from the config file. This is how often the
        # class update function is called. The all_hosts argument gives the full list of hosts at
        # the start of the simulation which may be required to initialise sizes. The argument
        # all_cells will provide list of cells if simulation type is "RASTER"

        # Set update frequency
        self.update_freq = update_freq

        # Set intervention type. If this is set to "CONTINUOUS" then intervention adds a rate
        # structure and events happen at those rates. Otherwise any name can be used for
        # interventions at particular times.
        self.type = "CONTINUOUS"

        # Set size of required rate structure. Only needed if self.type is "CONTINUOUS"
        self.rate_size = len(all_hosts)

        # Any other intervention set up code required

    # This function must be implemented if the intervention type is "CONTINUOUS". This defines the
    # events that must be carried out once a continuous rate has been selected from the rate
    # structure. Options for event_type are "CULL"
    def action(self, all_hosts, time, event_id, all_cells=None):
        # Function must return list of tuples: [(host_id, event_type)]
        event_list = []
        return event_list

    # This function is called every update_freq in time, and also after every event if
    # UpdateOnAllEvents is True. For "CONTINUOUS" interventions this must return a list of rate
    # changes to insert into rate structure: [(rate_id, new_rate)]. For any other intervention
    # return events to carry out: [(host_id, event_type)]. The argument after_event will be event
    # passed from event handler if UpdateOnAllEvents is True and this function is being called after
    # an event (rather than at a specified update_freq time). Options for event_type are "CULL".
    # Argument get_rate_fn takes position in rate structure and returns the rate currently stored
    # there if type is "CONTINUOUS"
    def update(self, all_hosts, time, all_cells=None, after_event=None, get_rate_fn=None,
               initial=False):
        if self.type == "CONTINUOUS":
            rate_list = []
            return rate_list
        else:
            event_list = []
            return event_list

    # Carry out any necessary operations at the end of an individual simulation. This function is
    # useful for storing information ready for finalised log file entry.
    def finalise(self):
        pass

    # Function to get log string to include in log file for set of simulations. Called after all
    # simulations have been carried out and finalise has been called for each.
    def output(self):
        log_str = ""
        return log_str
