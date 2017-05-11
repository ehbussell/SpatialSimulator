import config
import argparse
import eventhandling
import outputdata
import copy
import time as time_mod
import numpy as np
from host import Host
from ratesum import RateSum


def kernel(dist):
    if dist > 0:
        return 0.01/dist
    else:
        return 0


def calc_dist_kernel(hosts):
    nhosts = len(hosts)
    distances = np.zeros((nhosts, nhosts))
    kernel_vals = np.zeros((nhosts, nhosts))

    for i in range(nhosts):
        for j in range(i):
            dist = np.linalg.norm([hosts[i].x-hosts[j].x, hosts[i].y-hosts[j].y])
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

    time1 = time_mod.time()

    # Read in hosts
    params['init_hosts'] = config.read_hosts()
    params['nhosts'] = len(params['init_hosts'])

    # Setup
    params['next_state'] = setup(params)
    params['kernel'] = kernel
    params['init_region_summary'] = [{key: 0 for key in list(params['Model'])}
                                     for _ in range(params['NRegions'])]

    if params['RateStructure-Infection'] == "ratesum":
        inf_rates = RateSum(params['nhosts'])
    else:
        raise ValueError("Invalid rate structure - infection events!")

    if params['RateStructure-Advance'] == "ratesum":
        adv_rates = RateSum(params['nhosts'])
    else:
        raise ValueError("Invalid rate structure - advance events!")

    all_rates = [inf_rates, adv_rates]

    # Setup initial rates
    if params['CacheKernel'] is True:
        distances, kernel_vals = calc_dist_kernel(params['init_hosts'])
        params['kernel_vals'] = kernel_vals
        params['distances'] = distances
        do_event = eventhandling.do_event_cached

    else:
        do_event = eventhandling.do_event_uncached

    params['init_inf_rates'] = np.zeros(params['nhosts'])
    params['init_adv_rates'] = np.zeros(params['nhosts'])
    for i in range(params['nhosts']):
        current_state = params['init_hosts'][i].state
        params['init_region_summary'][params['init_hosts'][i].reg][current_state] += 1
        if current_state in "ECDI":
            params['init_adv_rates'][i] = params[current_state + 'AdvRate']
            if current_state in "CI":
                for j in range(params['nhosts']):
                    if params['init_hosts'][j].state == "S":
                        if params['CacheKernel'] is True:
                            params['init_inf_rates'][j] += params['kernel_vals'][j, i]
                        else:
                            params['init_inf_rates'][j] += kernel(np.linalg.norm([
                                params['init_hosts'][i].x-params['init_hosts'][j].x,
                                params['init_hosts'][i].y-params['init_hosts'][j].y]))

    print("Initial setup complete")

    for iteration in range(params['NIterations']):

        run_params = {}
        run_params['all_events'] = []
        run_params['region_summary'] = copy.deepcopy(params['init_region_summary'])
        run_params['summary_dump'] = []

        all_hosts = copy.deepcopy(params['init_hosts'])

        inf_rates.zero_rates()
        adv_rates.zero_rates()

        for i in range(params['nhosts']):
            inf_rates.insert_rate(i, params['init_inf_rates'][i])
            adv_rates.insert_rate(i, params['init_adv_rates'][i])

        time2 = time_mod.time()

        time = 0
        run_params['summary_dump'].append((time, copy.deepcopy(params['init_region_summary'])))
        nextSummaryDumpTime = time + params['SummaryOutputFreq']

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

            while time >= nextSummaryDumpTime:
                run_params['summary_dump'].append((nextSummaryDumpTime,
                                                  copy.deepcopy(run_params['region_summary'])))
                nextSummaryDumpTime += params['SummaryOutputFreq']

            if time > params['FinalTime']:
                break

            if select_rate < inf_rates.get_total_rate():
                eventID = inf_rates.select_event(select_rate)
                do_event(eventID, all_hosts, all_rates, params, run_params, time)
            elif select_rate < adv_rates.get_total_rate():
                eventID = adv_rates.select_event(select_rate - inf_rates.get_total_rate())
                do_event(eventID, all_hosts, all_rates, params, run_params, time)

        run_params['summary_dump'].append((nextSummaryDumpTime,
                                          copy.deepcopy(run_params['region_summary'])))

        time3 = time_mod.time()

        print(time2-time1, time3-time2, time3-time1)

        outputdata.output_data_hosts(all_hosts, params, iteration=iteration)
        outputdata.output_data_events(run_params, iteration=iteration)
        outputdata.output_data_summary(params, run_params, iteration=iteration)


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
