"""Tools for generating interevention scripts."""

import pprint
from IndividualSimulator.code.interventions import ContRegionRemoval


def create_continuous_region_removal_script(
        filename, nregions=3, priorities=[2, 1, 0], budget=200, removal_rate=1.0):
    """Create a region prioritised, continuous removal interevention, with budget constraint."""

    config_params = {
        'nregions': nregions,
        'priorities': priorities,
        'budget': budget,
        'removal_rate': removal_rate,
    }

    with open(ContRegionRemoval.__file__, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith("config_params"):
            lines[i] = "config_params = " + pprint.pformat(config_params, compact=False) + "\n"

    with open(filename, "w") as f:
        f.writelines(lines)
