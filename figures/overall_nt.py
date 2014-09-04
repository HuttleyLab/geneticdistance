from __future__ import division

from numpy import array, std, sqrt
import sys
import os

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r as R

import lib

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.2-dev'

def print_mean_and_ci(r, x):
    print r['name'] + ('+Gamma' if r['with_rate'] else '')
    print 'Samples', len(x)
    x = array(x)
    m, s = x.mean(), std(x)/sqrt(len(x))
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

