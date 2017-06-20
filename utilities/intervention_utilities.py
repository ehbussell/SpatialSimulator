"""Tools for generating interevention scripts."""

import pprint


def create_continuous_region_removal_script(
        filename, nregions=3, priorities=[2, 1, 0], budget=200, removal_rate=1.0):
    """Create a region prioritised, continuous removal interevention, with budget constraint."""

    config_params = {
        'nregions': nregions,
        'priorities': priorities,
        'budget': budget,
        'removal_rate': removal_rate,
    }

    with open(filename, "w") as f:
        f.write(
            "from IndividualSimulator.code.interventions.ContRegionRemoval import Intervention")
        f.write("\n\n")
        f.write("config_params = ")
        f.write(pprint.pformat(config_params, compact=False))
