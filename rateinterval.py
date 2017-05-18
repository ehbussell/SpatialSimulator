import numpy as np


class RateInterval:

    def __init__(self, size):
        self.interval_length = int(np.sqrt(size))
        if self.interval_length < 1:
            self.interval_length = 1

        self.n_intervals = int(np.ceil(size / self.interval_length))

        self.padded_length = self.interval_length * self.n_intervals

        self.zero_rates()

    def insert_rate(self, pos, rate):
        rate_change = rate - self.sub_individ_rates[pos]

        self.sub_individ_rates[pos] = rate

        interval_num = int(pos/self.interval_length)
        if interval_num < self.flag_changes:
            self.flag_changes = interval_num

        self.super_individ_rates[interval_num] += rate_change

        self.totrate += rate_change

    def get_rate(self, pos):
        return self.sub_individ_rates[pos]

    def select_event(self, rate):
        intervalID = self._interval_search_super(rate)

        rate_marginal = rate
        if intervalID != 0:
            rate_marginal = rate - self.super_sum_rates[intervalID - 1]

        return self._interval_search_sub(rate_marginal, intervalID)

    def get_total_rate(self):
        if self.flag_changes < self.n_intervals:
            self._sum_super_rates()
            self.totrate = self.super_sum_rates[self.n_intervals - 1]

        return self.totrate

    def full_resum(self):
        pass

    def zero_rates(self):
        self.sub_individ_rates = np.zeros(self.padded_length)
        self.sub_sum_rates = np.zeros(self.padded_length)

        self.super_individ_rates = np.zeros(self.n_intervals)
        self.super_sum_rates = np.zeros(self.n_intervals)

        self.flag_changes = self.n_intervals

        self.totrate = 0.0

    def _sum_super_rates(self):
        val = self.super_individ_rates[0]
        self.super_sum_rates[0] = val

        for i in range(1, self.n_intervals):
            val += self.super_individ_rates[i]
            self.super_sum_rates[i] = val

    def _sum_interval_rate(self, intervalID):
        first_idx = intervalID*self.interval_length
        last_idx = first_idx + self.interval_length - 1
        val = self.sub_individ_rates[first_idx]
        self.sub_sum_rates[first_idx] = val

        for i in range(1, self.interval_length):
            val += self.sub_individ_rates[first_idx + i]
            self.sub_sum_rates[first_idx + i] = val

        self.super_individ_rates[intervalID] = self.sub_sum_rates[last_idx]

    def _interval_search_sub(self, rate, intervalID):
        self._sum_interval_rate(intervalID)

        nlow = intervalID*self.interval_length
        nhigh = nlow + self.interval_length - 1

        while (nhigh - nlow) > 1:
            nmid = int(nlow + (nhigh-nlow)/2)
            if self.sub_sum_rates[nmid] > rate:
                nhigh = nmid
            else:
                nlow = nmid

        if self.sub_sum_rates[nlow] > rate:
            return nlow
        else:
            return nhigh

    def _interval_search_super(self, rate):
        nlow = 0
        nhigh = self.n_intervals - 1

        while (nhigh - nlow) > 1:
            nmid = int(nlow + (nhigh-nlow)/2)
            if self.super_sum_rates[nmid] > rate:
                nhigh = nmid
            else:
                nlow = nmid

        if self.super_sum_rates[nlow] > rate:
            return nlow
        else:
            return nhigh
