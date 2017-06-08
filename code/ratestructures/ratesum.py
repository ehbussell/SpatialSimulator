import numpy as np


class RateSum:

    def __init__(self, size):
        self.rates = np.zeros(size)
        self.totrate = np.sum(self.rates)
        self.nevents = size

    def insert_rate(self, pos, rate):
        if rate < 0:
            rate = 0
        rate_change = rate - self.rates[pos]
        self.rates[pos] = rate
        self.totrate += rate_change

    def get_rate(self, pos):
        return self.rates[pos]

    def select_event(self, rate):
        eventID = 0
        cum_rate = self.rates[0]

        while cum_rate < rate and eventID < self.nevents:
            eventID += 1
            cum_rate += self.rates[eventID]

        return eventID

    def get_total_rate(self):
        return self.totrate

    def full_resum(self):
        self.totrate = np.sum(self.rates)

    def zero_rates(self):
        self.rates = np.zeros(self.nevents)
        self.full_resum()
