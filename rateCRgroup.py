import numpy as np
from namedlist import namedlist

IndexStorage = namedlist("IndexStorage", "iIndex", default=None)
EventStorage = namedlist("EventStorage", "iEvent Rate", default=None)


class RateCRGroup:

    def __init__(self, parentCR, size, index, min_rate, max_rate):
        self.parentCR = parentCR
        self.n_events_max = size
        self.i_group_index = index
        self.min_rate = min_rate
        self.max_rate = max_rate

        self.location_to_index_map = [IndexStorage() for _ in range(self.n_events_max)]
        self.events = [EventStorage() for _ in range(self.n_events_max)]

        self.zero_rates()

    def insert_rate(self, pos, rate):
        eventIdx = self.location_to_index_map[pos].iIndex
        rate_change = rate - self.events[eventIdx].Rate

        self.events[eventIdx].Rate = rate
        self.totrate += rate_change

        self.parentCR._sub_group_report_rate(self.i_group_index, self.totrate)

    def get_rate(self, pos):
        return self.events[self.location_to_index_map[pos].iIndex].Rate

    def select_event(self, rate):
        selected_event = -1

        while selected_event < 0:
            random_rate = self.max_rate*np.random.random_sample()
            random_index = np.random.randint(0, self.n_events_active)

            if random_rate < self.events[random_index].Rate:
                selected_event = self.events[random_index].iEvent

        return selected_event

    def get_total_rate(self):
        return self.totrate

    def full_resum(self):
        pass

    def zero_rates(self):
        self.totrate = 0.0
        self.n_events_active = 0

        self.parentCR._sub_group_report_rate(self.i_group_index, self.totrate)

    def fill_zero_rates(self, n_rates_to_fill):
        self.totrate = 0.0
        self.n_events_active = n_rates_to_fill

        for i in range(n_rates_to_fill):
            self.events[i].iEvent = i
            self.events[i].Rate = 0.0
            self.location_to_index_map[i].iIndex = i

        self.parentCR._sub_group_report_rate(self.i_group_index, self.totrate)

    def insert_element_to_group(self, pos, rate):
        self.events[self.n_events_active].iEvent = pos
        self.events[self.n_events_active].Rate = rate

        self.location_to_index_map[pos].iIndex = self.n_events_active

        self.n_events_active += 1

        self.totrate += rate

        self.parentCR._sub_group_report_rate(self.i_group_index, self.totrate)

    def remove_element_from_group(self, pos):
        removal_index = self.location_to_index_map[pos].iIndex

        removed_rate = self.events[removal_index].Rate
        prev_last_pos = self.events[self.n_events_active-1].iEvent

        self.location_to_index_map[prev_last_pos].iIndex = removal_index

        self.events[removal_index].iEvent = prev_last_pos
        self.events[removal_index].Rate = self.events[self.n_events_active-1].Rate

        self.n_events_active -= 1
        self.totrate -= removed_rate

        self.parentCR._sub_group_report_rate(self.i_group_index, self.totrate)
