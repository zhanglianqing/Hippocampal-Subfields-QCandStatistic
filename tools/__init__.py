from .FSsubfieldsQC import post_SF_segment
from .utils.statistics import merge_table, Wrap_Statistic
tasks = {
    'QA': post_SF_segment,
    'merge':merge_table,
    'statistic':Wrap_Statistic
}


def create_task(tname):
    module = tasks.get(tname, None)
    return module