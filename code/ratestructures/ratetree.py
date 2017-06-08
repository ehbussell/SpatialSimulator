import numpy as np


class RateTree:

    def __init__(self, size):
        self.n_tree_levels = 1
        self.padded_length = 1
        while self.padded_length < size:
            self.padded_length *= 2
            self.n_tree_levels += 1

        self.rates = np.zeros(2*self.padded_length - 1)
        self.tree_levels = [None for _ in range(self.n_tree_levels + 2)]

        self._calculate_tree_levels()

        self.zero_rates()

    def insert_rate(self, pos, rate):
        level = 1
        loc = pos

        rate_change = rate - self.rates[self.tree_levels[level] + loc]

        self.totrate += rate_change

        while self.tree_levels[level] is not None:
            self.rates[self.tree_levels[level] + loc] += rate_change
            level += 1
            loc = loc >> 1

    def get_rate(self, pos):
        return self.rates[pos]

    def select_event(self, rate):
        level = self.n_tree_levels - 1

        idx = 0

        while self.tree_levels[level] is not None:
            idx *= 2

            leftRate = self.rates[self.tree_levels[level] + idx]
            if leftRate <= rate:
                idx += 1
                rate -= leftRate

            level -= 1

        return idx

    def get_total_rate(self):
        return self.totrate

    def full_resum(self):
        level_length = self.padded_length/2

        for i in range(2, self.n_tree_levels + 1):
            for j in range(int(level_length)):
                self.rates[self.tree_levels[i] + j] = self.rates[self.tree_levels[i-1] + 2*j] + \
                    self.rates[self.tree_levels[i-1] + 2*j + 1]
            level_length /= 2

    def zero_rates(self):
        for i in range(2*self.padded_length - 1):
            self.rates[i] = 0

        self.totrate = 0

    def _calculate_tree_levels(self):
        self.tree_levels[0] = None
        self.tree_levels[self.n_tree_levels + 1] = None
        level_length = self.padded_length
        level_start = 0
        for i in range(1, self.n_tree_levels+1):
            self.tree_levels[i] = int(level_start)
            level_start += level_length
            level_length /= 2
