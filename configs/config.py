import os
import re
import argparse
from datetime import datetime
from .yacs import CfgNode as CN


# -----------------------------------------------------------------------------
# global setting
# -----------------------------------------------------------------------------
cfg = CN()
cfg.task_total = ['QA','merge','statistic']
#
#

# -----------------------------------------------------------------------------
# task setting
# -----------------------------------------------------------------------------
cfg.task = CN()


# -----------------------------------------------------------------------------
# other setting
# -----------------------------------------------------------------------------
cfg.other = CN()


def make_cfg_args():
    parser = argparse.ArgumentParser(description='Generate QA snapshots and summary tables')
    parser.add_argument('--data_root', type=str, metavar='PATH', help='path to $SUBJECTS_DIR')
    parser.add_argument('--task_name', default='QA', type=str, help='which tasks?')
    parser.add_argument('--input_table',type=str, metavar='PATH', help='path to clinical_data.csv')
    parser.add_argument('--SF_volume',type=str, metavar='PATH', help='path to SF_volume.csv')
    parser.add_argument('--output',type=str, metavar='PATH', help='path to output directory')
    parser.add_argument('--design',type=str, metavar='PATH', help='path to design.csv')
    parser.add_argument('--Rcmd',default='Rscript',type=str,help='commend to run Rscript. e.g. C:/PROGRA~1/R/R-41~1.1/bin/Rscript')
    parser.add_argument('--n_cpus',default='all',type=str,help='number of cores for processing')
    parser.add_argument('--no_plot',action="store_false",help="don't generate qa plot,just generate volume table")
    args = parser.parse_args()
    assert args.task_name in cfg.task_total, 'error task name ({})'.format(args.task_name)
    import multiprocessing
    if args.n_cpus == 'all':
        args.n_cpus = multiprocessing.cpu_count()
    else:
        args.n_cpus = int(args.n_cpus)

    return cfg.clone(), args