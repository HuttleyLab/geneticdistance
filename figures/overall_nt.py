from __future__ import division

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r as R

import numpy as np
from scipy.stats import sem
import sys
import os

def print_mean_and_ci(r, x):
    print r['name'] + ('+Gamma' if r['with_rate'] else '')
    print 'Samples', len(x)
    x = np.array(x)
    m, s = x.mean(), sem(x)
    print 'mean', 'stdev', 'ci'
    print m, s, m-1.96*s, m+1.96*s


def plot_ENS_hexbin(stats, pair=None, output_file=None, **kw):
    pair = list(pair)
    x = [(r[0]['EN'][pair[0]]/r[0]['EN'][pair[1]])
            for t in stats for r in stats[t]]
    y = [(r[1]['EN'][pair[0]]/r[1]['EN'][pair[1]]) 
            for t in stats for r in stats[t] if r[1]['EN'][pair[1]] > 1e-9]
    z = [(r[3]['EN'][pair[0]]/r[3]['EN'][pair[1]]) 
            for t in stats for r in stats[t]]


    print pair[0] + ' / ' + pair[1]
    print_mean_and_ci(r[0], x)
    print_mean_and_ci(r[1], y)
    print_mean_and_ci(r[2], z)

def main():
    from results import get_results
    stats, args = get_results()

    plot_ENS_hexbin(stats, **vars(args)) # vb_two
        
if __name__ == '__main__':
    main()

