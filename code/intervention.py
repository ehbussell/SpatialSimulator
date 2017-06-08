class Intervention:
    def __init__(self, update_freq, parent_sim):
        """Setup any necessities for interventions."""
        self.update_freq = update_freq
        self.parent_sim = parent_sim

    def update(self):
        """Carry out interventions."""
        pass


class ContinuousIntervention:
    def __init__(self, update_freq, parent_sim, rate_handler):
        """Setup any necessities for interventions."""
        self.update_freq = update_freq
        self.parent_sim = parent_sim
        self.rate_handler = rate_handler

    def update(self):
        """Adjust rates and factors."""
        pass

    def action(self, eventID):
        """Carry out event."""
        pass
