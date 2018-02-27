import unittest
import os
import numpy as np
from collections import Counter
import scipy.stats
import matplotlib.pyplot as plt
from IndividualSimulator.code.ratestructures.ratesum import RateSum
from IndividualSimulator.code.ratestructures.ratetree import RateTree
from IndividualSimulator.code.ratestructures.rateinterval import RateInterval
from IndividualSimulator.code.ratestructures.rateCR import RateCR

def initialise_rates(rate_struct, size, rates=None):
    """Initialise rate structure rates, by default random numbers."""

    if rates is None:
        rates = np.random.rand(size)

    for i, rate in enumerate(rates):
        rate_struct.insert_rate(i, rate)

    return rates

def run_selections(rate_struct, size, niters):
    """Run niters rate selections from rate_struct, returning number of selections for each event"""

    total_rate = rate_struct.get_total_rate()

    cdict = Counter()
    for i in range(niters):
        select_rate = np.random.rand()*total_rate
        selected_id = rate_struct.select_event(select_rate)
        cdict[selected_id] += 1

    all_n_selected = [cdict[i] for i in range(size)]

    return all_n_selected

class RateSumTests(unittest.TestCase):
    """Test that rate sum structure performs correctly."""

    def setUp(self):
        self.size = 1000
        self.rate_struct = RateSum(self.size)

    def test_get_insert(self):
        "Test RateSum structure get/insert rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        for i, rate in enumerate(new_rates):
            self.assertEqual(self.rate_struct.get_rate(i), rate)

    def test_total_rate(self):
        "Test RateSum structure get total rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.assertAlmostEqual(self.rate_struct.get_total_rate(), np.sum(new_rates))

    def test_zero_rates(self):
        "Test RateSum structure zero rates functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.rate_struct.zero_rates()

        for i in range(self.size):
            self.assertEqual(self.rate_struct.get_rate(i), 0.0)

        self.assertEqual(self.rate_struct.get_total_rate(), 0.0)

    def test_select_rate(self):
        "Test RateSum structure select rate functions."""

        # Check uniform rates
        new_rates = initialise_rates(self.rate_struct, self.size, np.ones(self.size))
        niters = int(1e5)
        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        mu = niters / self.size
        x = np.arange(scipy.stats.poisson.ppf(0.001, mu), scipy.stats.poisson.ppf(0.999, mu))

        plt.style.use("ggplot")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(all_n_selected, normed=True, bins=30, color="blue", alpha=0.3)
        ax.plot(x, scipy.stats.poisson.pmf(x, mu), 'r--', lw=2, label='poisson pmf')
        ax.set_xlabel("Number of Rate Selections")
        ax.set_ylabel("Frequency")
        ax.set_title("RateSum Test Results")
        fig.savefig(os.path.join("testing", "RateSumTest.png"))

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)

        # Check random rates
        new_rates = initialise_rates(self.rate_struct, self.size)

        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)


class RateTreeTests(unittest.TestCase):
    """Test that rate tree structure performs correctly."""

    def setUp(self):
        self.size = 1000
        self.rate_struct = RateTree(self.size)

    def test_get_insert(self):
        "Test RateTree structure get/insert rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        for i, rate in enumerate(new_rates):
            self.assertEqual(self.rate_struct.get_rate(i), rate)

    def test_total_rate(self):
        "Test RateTree structure get total rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.assertAlmostEqual(self.rate_struct.get_total_rate(), np.sum(new_rates))

    def test_zero_rates(self):
        "Test RateTree structure zero rates functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.rate_struct.zero_rates()

        for i in range(self.size):
            self.assertEqual(self.rate_struct.get_rate(i), 0.0)

        self.assertEqual(self.rate_struct.get_total_rate(), 0.0)

    def test_select_rate(self):
        "Test RateTree structure select rate functions."""

        # Check uniform rates
        new_rates = initialise_rates(self.rate_struct, self.size, np.ones(self.size))
        niters = int(1e5)
        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        mu = niters / self.size
        x = np.arange(scipy.stats.poisson.ppf(0.001, mu), scipy.stats.poisson.ppf(0.999, mu))

        plt.style.use("ggplot")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(all_n_selected, normed=True, bins=30, color="blue", alpha=0.3)
        ax.plot(x, scipy.stats.poisson.pmf(x, mu), 'r--', lw=2, label='poisson pmf')
        ax.set_xlabel("Number of Rate Selections")
        ax.set_ylabel("Frequency")
        ax.set_title("RateTree Test Results")
        fig.savefig(os.path.join("testing", "RateTreeTest.png"))

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)

        # Check random rates
        new_rates = initialise_rates(self.rate_struct, self.size)

        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)


class RateIntervalTests(unittest.TestCase):
    """Test that rate interval structure performs correctly."""

    def setUp(self):
        self.size = 1000
        self.rate_struct = RateInterval(self.size)

    def test_get_insert(self):
        "Test RateInterval structure get/insert rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        for i, rate in enumerate(new_rates):
            self.assertEqual(self.rate_struct.get_rate(i), rate)

    def test_total_rate(self):
        "Test RateInterval structure get total rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.assertAlmostEqual(self.rate_struct.get_total_rate(), np.sum(new_rates))

    def test_zero_rates(self):
        "Test RateInterval structure zero rates functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.rate_struct.zero_rates()

        for i in range(self.size):
            self.assertEqual(self.rate_struct.get_rate(i), 0.0)

        self.assertEqual(self.rate_struct.get_total_rate(), 0.0)

    def test_select_rate(self):
        "Test RateInterval structure select rate functions."""

        # Check uniform rates
        new_rates = initialise_rates(self.rate_struct, self.size, np.ones(self.size))
        niters = int(1e5)
        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        mu = niters / self.size
        x = np.arange(scipy.stats.poisson.ppf(0.001, mu), scipy.stats.poisson.ppf(0.999, mu))

        plt.style.use("ggplot")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(all_n_selected, normed=True, bins=30, color="blue", alpha=0.3)
        ax.plot(x, scipy.stats.poisson.pmf(x, mu), 'r--', lw=2, label='poisson pmf')
        ax.set_xlabel("Number of Rate Selections")
        ax.set_ylabel("Frequency")
        ax.set_title("RateInterval Test Results")
        fig.savefig(os.path.join("testing", "RateIntervalTest.png"))

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)

        # Check random rates
        new_rates = initialise_rates(self.rate_struct, self.size)

        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)


class RateCRTests(unittest.TestCase):
    """Test that rate CR structure performs correctly."""

    def setUp(self):
        self.size = 1000
        self.rate_struct = RateTree(self.size)

    def test_get_insert(self):
        "Test RateCR structure get/insert rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        for i, rate in enumerate(new_rates):
            self.assertEqual(self.rate_struct.get_rate(i), rate)

    def test_total_rate(self):
        "Test RateCR structure get total rate functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.assertAlmostEqual(self.rate_struct.get_total_rate(), np.sum(new_rates))

    def test_zero_rates(self):
        "Test RateCR structure zero rates functions."""

        new_rates = initialise_rates(self.rate_struct, self.size)

        self.rate_struct.zero_rates()

        for i in range(self.size):
            self.assertEqual(self.rate_struct.get_rate(i), 0.0)

        self.assertEqual(self.rate_struct.get_total_rate(), 0.0)

    def test_select_rate(self):
        "Test RateCR structure select rate functions."""

        # Check uniform rates
        new_rates = initialise_rates(self.rate_struct, self.size, np.ones(self.size))
        niters = int(1e5)
        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        mu = niters / self.size
        x = np.arange(scipy.stats.poisson.ppf(0.001, mu), scipy.stats.poisson.ppf(0.999, mu))

        plt.style.use("ggplot")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(all_n_selected, normed=True, bins=30, color="blue", alpha=0.3)
        ax.plot(x, scipy.stats.poisson.pmf(x, mu), 'r--', lw=2, label='poisson pmf')
        ax.set_xlabel("Number of Rate Selections")
        ax.set_ylabel("Frequency")
        ax.set_title("RateCR Test Results")
        fig.savefig(os.path.join("testing", "RateCRTest.png"))

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)

        # Check random rates
        new_rates = initialise_rates(self.rate_struct, self.size)

        all_n_selected = run_selections(self.rate_struct, self.size, niters)

        f_obs = [all_n_selected[i]/niters for i in range(self.size)]
        f_exp = [self.rate_struct.get_rate(i)/self.rate_struct.get_total_rate()
                 for i in range(self.size)]

        chi, pval = scipy.stats.chisquare(f_obs, f_exp)

        self.assertTrue(pval > 0.1)
