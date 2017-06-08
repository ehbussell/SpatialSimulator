import numpy as np
from .ratestructures.ratesum import RateSum
from .ratestructures.rateinterval import RateInterval
from .ratestructures.ratetree import RateTree
from .ratestructures.rateCR import RateCR


class RateHandler:

    def __init__(self, parent_sim):
        self.parent_sim = parent_sim
        self.params = self.parent_sim.params

        if self.params['RateStructure-Infection'] == "ratesum":
            self.inf_rates = RateSum(self.params['nhosts'])
        elif self.params['RateStructure-Infection'] == "rateinterval":
            self.inf_rates = RateInterval(self.params['nhosts'])
        elif self.params['RateStructure-Infection'] == "ratetree":
            self.inf_rates = RateTree(self.params['nhosts'])
        elif self.params['RateStructure-Infection'] == "rateCR":
            self.inf_rates = RateCR(self.params['nhosts'], 0.125,
                                    self.params['nhosts']*self.params['nhosts'])
        else:
            raise ValueError("Invalid rate structure - infection events!")

        if self.params['RateStructure-Advance'] == "ratesum":
            self.adv_rates = RateSum(self.params['nhosts'])
        elif self.params['RateStructure-Advance'] == "rateinterval":
            self.adv_rates = RateInterval(self.params['nhosts'])
        elif self.params['RateStructure-Advance'] == "ratetree":
            self.adv_rates = RateTree(self.params['nhosts'])
        elif self.params['RateStructure-Advance'] == "rateCR":
            self.adv_rates = RateCR(self.params['nhosts'], 0.125,
                                    self.params['nhosts']*self.params['nhosts'])
        else:
            raise ValueError("Invalid rate structure - advance events!")

        self.event_types = ["Infection", "Advance"]
        self.all_rates = {"Infection": self.inf_rates, "Advance": self.adv_rates}

        self.n_types = range(len(self.event_types))

    def add_rate_struct(self, event_type, size, structure="ratesum"):
        if event_type in self.all_rates:
            raise ValueError("Event type already exists!")

        self.event_types.append(event_type)
        rate_structure = RateSum(size)
        rate_structure.zero_rates()
        self.all_rates[event_type] = rate_structure

        self.parent_sim.rate_factor.append(1)

    def get_total_rate(self):
        total_rate = np.sum([self.parent_sim.rate_factor[i] *
                             self.all_rates[rate_type].get_total_rate()
                             for i, rate_type in enumerate(self.event_types)])
        return total_rate

    def get_next_event(self):
        total_rates = [self.parent_sim.rate_factor[i] *
                       self.all_rates[rate_type].get_total_rate()
                       for i, rate_type in enumerate(self.event_types)]
        total_rate = np.sum(total_rates)
        if total_rate < 10e-10:
            return (total_rate, None, None)

        select_rate = np.random.random_sample()*total_rate

        cumulative_rate = 0.0

        for i in range(len(self.event_types)):
            group_rate = total_rates[i]
            if select_rate < cumulative_rate + group_rate:
                return (total_rate, self.event_types[i],
                        self.all_rates[self.event_types[i]].select_event(
                            (select_rate - cumulative_rate)/self.parent_sim.rate_factor[i]))
            else:
                cumulative_rate += group_rate

        return (total_rate, None, None)

    def zero_rates(self):
        for event_type in self.event_types:
            self.all_rates[event_type].zero_rates()

    def insert_rate(self, hostID, rate, rate_type):
        return self.all_rates[rate_type].insert_rate(hostID, rate)

    def get_rate(self, hostID, rate_type):
        return self.all_rates[rate_type].get_rate(hostID)
