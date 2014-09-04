from __future__ import division

import numpy as np
from collections import defaultdict, Counter
from results import gs_p, get_results
import sys

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, IntVector, StrVector, globalenv
from rpy2.robjects import r as R

import lib

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def filter_akaike(rs, model):
    ak = {r['name']+('+Gamma' if r['with_rate'] else ''):-2*r['ll']+2*r['df']
            for r in rs}
    other_model = 'GTR' if model == 'GTR+Gamma' else 'GTR+Gamma'
    return ak[model] < ak[other_model]

def plot_bar(stats, pair=None, model=None, **kw):
    names = [r['name'] for r in stats.values()[0][0]]
    with_rates = [r['with_rate'] for r in stats.values()[0][0]]
    names = [n+('+Gamma' if w else '') for n, w in zip(names, with_rates)]

    by_dir = defaultdict(list)
    for triad in stats:
        for r in stats[triad]:
            if model is None or filter_akaike(r, model):
                by_dir[r[0]['from_directory']].append(r)

    for d in by_dir:
        by_dir[d] = zip(*[[sum(_r['EN'][p] for p in pair) for _r in r] 
            for r in by_dir[d]])

    runs = []
    g_stats = []
    data = []
    alpha = 0
    for d, v in by_dir.items():
        if 'exons' in d.split('/'):
            dataset = 'Nuclear'
        elif 'mtDNA' in d.split('/'):
            dataset = 'Mitochondrial'
        else:
            dataset = 'Microbial'
        print dataset
        for j, ens in enumerate(v):
            print names[j], np.mean(ens), np.std(ens)
        print 'Samples', len(ens)

    

def main():
    stats, args = get_results()

    plot_bar(stats, **vars(args))

    return 0

if __name__ == '__main__':
    sys.exit(main())
