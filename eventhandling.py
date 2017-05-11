import numpy as np


def do_event_cached(hostID, all_hosts, all_rates, params, run_params, time):
    """Carry out event on given host, and update rates."""
    inf_rates, adv_rates = all_rates

    old_state = all_hosts[hostID].state
    new_state = params['next_state'](old_state)
    all_hosts[hostID].state = new_state

    if old_state == "S":
        inf_rates.insert_rate(hostID, 0.0)
    if new_state in "ECDI":
        adv_rates.insert_rate(hostID, params[new_state + 'AdvRate'])
    if new_state == "R":
        adv_rates.insert_rate(hostID, 0.0)
        # Distribute rate changes
        for i in range(len(all_hosts)):
            if all_hosts[i].state == "S":
                old_rate = inf_rates.get_rate(i)
                new_rate = old_rate - params['kernel_vals'][i, hostID]
                inf_rates.insert_rate(i, new_rate)

    if new_state in "CI":
        # Distribute rate changes
        for i in range(len(all_hosts)):
            if all_hosts[i].state == "S":
                old_rate = inf_rates.get_rate(i)
                new_rate = old_rate + params['kernel_vals'][i, hostID]
                inf_rates.insert_rate(i, new_rate)

    all_hosts[hostID].jump_times[new_state] = time

    run_params['all_events'].append((time, hostID, old_state, new_state))
    run_params['region_summary'][all_hosts[hostID].reg][old_state] -= 1
    run_params['region_summary'][all_hosts[hostID].reg][new_state] += 1


def do_event_uncached(hostID, all_hosts, all_rates, params, run_params, time):
    """Carry out event on given host, and update rates."""
    inf_rates, adv_rates = all_rates

    old_state = all_hosts[hostID].state
    new_state = params['next_state'](old_state)
    all_hosts[hostID].state = new_state

    if old_state == "S":
        inf_rates.insert_rate(hostID, 0.0)
    if new_state in "ECDI":
        adv_rates.insert_rate(hostID, params[new_state + 'AdvRate'])
    if new_state == "R":
        adv_rates.insert_rate(hostID, 0.0)
        # Distribute rate changes
        for i in range(len(all_hosts)):
            if all_hosts[i].state == "S":
                old_rate = inf_rates.get_rate(i)
                new_rate = old_rate - params['kernel'](
                    np.linalg.norm([all_hosts[i].x-all_hosts[hostID].x,
                                    all_hosts[i].y-all_hosts[hostID].y]))
                inf_rates.insert_rate(i, new_rate)

    if new_state in "CI":
        # Distribute rate changes
        for i in range(len(all_hosts)):
            if all_hosts[i].state == "S":
                old_rate = inf_rates.get_rate(i)
                new_rate = old_rate + params['kernel'](
                    np.linalg.norm([all_hosts[i].x-all_hosts[hostID].x,
                                    all_hosts[i].y-all_hosts[hostID].y]))
                inf_rates.insert_rate(i, new_rate)

    all_hosts[hostID].jump_times[new_state] = time

    run_params['all_events'].append((time, hostID, old_state, new_state))
    run_params['region_summary'][all_hosts[hostID].reg][old_state] -= 1
    run_params['region_summary'][all_hosts[hostID].reg][new_state] += 1
