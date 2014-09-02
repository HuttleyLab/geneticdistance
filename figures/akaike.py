from __future__ import division

import lib

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, IntVector, StrVector, globalenv
from rpy2.robjects import r as R

import numpy as np

from collections import defaultdict, Counter
from results import gs_p, get_results
import sys

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def akaike_order(rs):
    ak = [(-2*r['ll']+2*r['df'],r['name']+('+Gamma' if r['with_rate'] else ''))
            for r in rs]
    ak.sort()
    return tuple(zip(*ak)[1])

def plot_bar(stats, **kw):
    names = [r['name'] for r in stats.values()[0][0]]
    with_rates = [r['with_rate'] for r in stats.values()[0][0]]
    names = [n+('+Gamma' if w else '') for n, w in zip(names, with_rates)]

    akaike = defaultdict(Counter)
    for triad in stats:
        for r in stats[triad]:
            akaike[r[0]['from_directory']][akaike_order(r)] += 1

    for d in akaike:
        if 'exons' in d.split('/'):
            dataset = 'Nuclear'
        elif 'mtDNA' in d.split('/'):
            dataset = 'Mitochondrial'
        else:
            dataset = 'Microbial'
        print dataset
        for order in akaike[d]:
            print order, akaike[d][order]

def main():
    stats, args = get_results()

    plot_bar(stats, **vars(args))

    return 0

if __name__ == '__main__':
    sys.exit(main())
