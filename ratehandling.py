import numpy as np
from ratesum import RateSum
from rateinterval import RateInterval
from ratetree import RateTree
from rateCR import RateCR


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

        self.all_rates = [self.inf_rates, self.adv_rates]
        self.event_types = ["Infection", "Advance"]

    def get_total_rate(self):
        total_rate = np.sum([self.parent_sim.rate_factor[i]*rate.get_total_rate()
                             for i, rate in enumerate(self.all_rates)])
        return total_rate

    def get_next_event(self):
        total_rates = [self.parent_sim.rate_factor[i]*rate.get_total_rate()
                       for i, rate in enumerate(self.all_rates)]
        total_rate = np.sum(total_rates)

        select_rate = np.random.random_sample()*total_rate

        cumulative_rate = 0.0

        for i in range(len(self.event_types)):
            group_rate = total_rates[i]
            if select_rate < cumulative_rate + group_rate:
                return (total_rate, self.event_types[i],
                        self.all_rates[i].select_event(select_rate - cumulative_rate))
            else:
                cumulative_rate += group_rate

        return (total_rate, None, None)

    def zero_rates(self):
        for rate in self.all_rates:
            rate.zero_rates()

    def insert_rate(self, hostID, rate, rate_type):
        if rate_type == "Infection":
            self.all_rates[0].insert_rate(hostID, rate)
        elif rate_type == "Advance":
            self.all_rates[1].insert_rate(hostID, rate)
        else:
            raise ValueError("Unrecognised event type")

    def get_rate(self, hostID, rate_type):
        if rate_type == "Infection":
            return self.all_rates[0].get_rate(hostID)
        elif rate_type == "Advance":
            return self.all_rates[1].get_rate(hostID)
        else:
            raise ValueError("Unrecognised event type")
