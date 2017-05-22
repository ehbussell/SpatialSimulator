import numpy as np
from namedlist import namedlist
from rateCRgroup import RateCRGroup
from ratetree import RateTree


IndexStorage = namedlist("IndexStorage", "iGroup", default=None)


class RateCR:

    def __init__(self, size, min_rate, max_rate):
        self.n_events_stored = size

        log_min_rate = np.log2(min_rate)
        self.i_offset_min_rate = int(np.floor(log_min_rate))
        self.min_rate = np.power(2.0, self.i_offset_min_rate)

        log_max_rate = np.log2(max_rate)
        self.i_offset_max_rate = int(np.ceil(log_max_rate))
        self.max_rate = np.power(2.0, self.i_offset_max_rate)

        self.n_groups_normal = self.i_offset_max_rate - self.i_offset_min_rate

        self.n_groups = self.n_groups_normal

        self.i_group_epsilon = self.n_groups
        self.n_groups += 1
        self.i_group_omega = self.n_groups
        self.n_groups += 1
        self.i_group_zero = self.n_groups
        self.n_groups += 1

        self.n_groups_special = self.n_groups - self.n_groups_normal

        self.group_rates = RateTree(self.n_groups)

        self.groups = []
        n_factor = 1
        for i in range(self.n_groups_normal):
            lower_rate = self.min_rate*n_factor
            self.groups.append(RateCRGroup(
                self, self.n_events_stored, i, lower_rate, lower_rate*2))
            n_factor *= 2

        self.groups.append(RateCRGroup(
            self, self.n_events_stored, self.i_group_epsilon, 0.0, self.min_rate))
        self.groups.append(RateCRGroup(
            self, self.n_events_stored, self.i_group_omega, self.max_rate, self.max_rate*32))
        self.groups.append(RateCRGroup(
            self, self.n_events_stored, self.i_group_zero, 0.0, 0.0))

        self.location_to_storage_map = [IndexStorage() for _ in range(self.n_events_stored)]

        self.zero_rates()

    def insert_rate(self, pos, rate):
        target_group = self._get_group_id_from_rate(rate)

        current_group = self.location_to_storage_map[pos].iGroup
        if current_group != target_group:
            self.groups[current_group].remove_element_from_group(pos)
            self.groups[target_group].insert_element_to_group(pos, rate)
            self.location_to_storage_map[pos].iGroup = target_group
        else:
            self.groups[target_group].insert_rate(pos, rate)

    def get_rate(self, pos):
        return self.groups[self.location_to_storage_map[pos].iGroup].get_rate(pos)

    def select_event(self, rate):
        cumulative_rate = 0.0

        for i in range(self.n_groups):
            group_rate = self.groups[i].get_total_rate()
            if rate < cumulative_rate + group_rate:
                return self.groups[i].select_event(rate)
            else:
                cumulative_rate += group_rate

    def get_total_rate(self):
        return self.group_rates.get_total_rate()

    def full_resum(self):
        pass

    def zero_rates(self):
        for i in range(self.n_groups):
            self.groups[i].zero_rates()

        self.groups[self.i_group_zero].fill_zero_rates(self.n_events_stored)

        self.group_rates.zero_rates()

        for i in range(self.n_events_stored):
            self.location_to_storage_map[i].iGroup = self.i_group_zero

    def _get_group_id_from_rate(self, rate):

        if rate > 0:
            if rate >= self.min_rate:
                if rate < self.max_rate:
                    log_rate = np.log2(rate)

                    return int(np.floor(log_rate)) - self.i_offset_min_rate
                else:
                    return self.i_group_omega
            else:
                return self.i_group_epsilon
        else:
            return self.i_group_zero

    def _sub_group_report_rate(self, group, rate):
        self.group_rates.insert_rate(group, rate)
