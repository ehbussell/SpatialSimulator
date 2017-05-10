import config
import argparse
import eventhandling
import numpy as np
from host import Host
from ratesum import RateSum


def kernel(dist):
    if dist > 0:
        return 1/dist
    else:
        return 0


def calc_dist_kernel(hosts):
    nhosts = len(hosts)
    distances = np.zeros((nhosts, nhosts))
    kernel_vals = np.zeros((nhosts, nhosts))

    for i in range(nhosts):
        for j in range(i):
            dist = np.linalg.norm([hosts[i].x-hosts[j].x,
                                   hosts[i].y-hosts[j].y])
            distances[i, j] = dist
            distances[j, i] = dist

            kval = kernel(dist)
            kernel_vals[i, j] = kval
            kernel_vals[j, i] = kval

    return (distances, kernel_vals)


def setup(params):
    states = list(params['Model'])

    def next_state_func(current_state):
        states_iter = iter(states)
        state = next(states_iter, None)
        while current_state != state:
            state = next(states_iter, None)
            if state is None:
                break
        return next(states_iter, None)

    return next_state_func


def run_epidemic(params):

    print(params)

    # Read in hosts
    all_hosts = config.read_hosts()
    nhosts = len(all_hosts)

    # Setup
    params['next_state'] = setup(params)

    if params['RateStructure-Infection'] == "ratesum":
        inf_rates = RateSum(nhosts)
    else:
        raise ValueError("Invalid rate structure - infection events!")

    if params['RateStructure-Advance'] == "ratesum":
        adv_rates = RateSum(nhosts)
    else:
        raise ValueError("Invalid rate structure - advance events!")

    all_rates = [inf_rates, adv_rates]

    # Setup initial rates
    if params['CacheKernel'] is True:
        distances, kernel_vals = calc_dist_kernel(all_hosts)
        params['kernel_vals'] = kernel_vals
        params['distances'] = distances
        do_event = eventhandling.do_event_cached

        inf_rates_tmp = np.zeros(nhosts)
        adv_rates_tmp = np.zeros(nhosts)
        for i in range(nhosts):
            current_state = all_hosts[i].state
            if current_state in "ECDI":
                adv_rates_tmp[i] = params[current_state + 'AdvRate']
                if current_state in "CI":
                    for j in range(nhosts):
                        if all_hosts[j].state == "S":
                            inf_rates_tmp[j] += params['kernel_vals'][j, i]

    else:
        do_event = eventhandling.do_event_uncached

        inf_rates_tmp = np.zeros(nhosts)
        adv_rates_tmp = np.zeros(nhosts)
        for i in range(nhosts):
            current_state = all_hosts[i].state
            if current_state in "ECDI":
                adv_rates_tmp[i] = params[current_state + 'AdvRate']
                if current_state in "CI":
                    for j in range(nhosts):
                        if all_hosts[j].state == "S":
                            inf_rates_tmp[j] += kernel(np.linalg.norm([
                                all_hosts[i].x-all_hosts[j].x,
                                all_hosts[i].y-all_hosts[j].y]))

    for i in range(nhosts):
        inf_rates.insert_rate(i, inf_rates_tmp[i])
        adv_rates.insert_rate(i, adv_rates_tmp[i])

    print("Initial setup complete")

    time = 0
    print_state(time, all_hosts)

    # Run gillespie loop
    while True:
        inf_rates.full_resum()  # As errors accumulate in total rate
        adv_rates.full_resum()  # As errors accumulate in total rate
        totRate = inf_rates.get_total_rate() + adv_rates.get_total_rate()
        if totRate <= 0:
            break

        nextTime = (-1.0/totRate)*np.log(np.random.random_sample())
        select_rate = np.random.random_sample()*totRate

        time += nextTime

        if time > params['FinalTime']:
            break

        if select_rate < inf_rates.get_total_rate():
            eventID = inf_rates.select_event(select_rate)
            do_event(eventID, all_hosts, all_rates, params)
        elif select_rate < adv_rates.get_total_rate():
            eventID = adv_rates.select_event(
                select_rate - inf_rates.get_total_rate())
            do_event(eventID, all_hosts, all_rates, params)

        print_state(time, all_hosts)


def print_state(time, all_hosts):
    print("Time: " + str(time))
    for host in all_hosts:
        print(host.state)
    print("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-c", "--configFile", default="config.ini",
                        help="Name of config file to use.")
    parser.add_argument("-k", "--keyFile", action="store_true",
                        help="Flag to generate keyfile.  "
                        "If present keyfile is created and program exits.")
    parser.add_argument("-d", "--defaultConfig", const="config.ini",
                        nargs="?", help="Flag to generate default config file."
                        "  If present file is created and program exits.  "
                        "Optionally specify filename.")
    args = parser.parse_args()

    if args.keyFile is True:
        config.write_keyfile()
        print("KeyFile generated.")

    if args.defaultConfig is not None:
        config.write_default_config(args.defaultConfig)
        print("Default config file generated.")

    if args.keyFile is False and args.defaultConfig is None:
        params = config.read_config_file(filename=args.configFile)
        run_epidemic(params)
