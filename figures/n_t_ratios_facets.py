from __future__ import division

import sys
import os

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv, StrVector, \
        FactorVector
from rpy2.robjects import r as R

import lib
from outliers import qcrop

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def label(row):
    if row['with_rate']:
        return row['name']+'+Gamma'
    return row['name']

def plot_n_t_ratios(stats, pair=None, output_file=None, **kw):
    pair = list(pair)
    x = []
    y = []
    facet = []
    labels = []
    for o in range(0, len(stats.values()[0][0]), 2):
        datalabel = label(stats.values()[0][0][o+1])
        labels.append(datalabel)
        samples = 0
        for triad in stats:
            for r in stats[triad]:
                x.append(r[o]['EN'][pair[0]]/r[o]['EN'][pair[1]])
                y.append(r[o+1]['EN'][pair[0]]/r[o+1]['EN'][pair[1]])
                facet.append(datalabel)
                samples += 1

    print 'Samples', samples

    xlab = pair[0] + ' over ' + pair[1] + ' Edge Length under General Model'
    ylab = pair[0] + ' over ' + pair[1] + \
            ' Edge Length under Stationary Models'

    globalenv['df'] = DataFrame({'x': FloatVector(x), 'y': FloatVector(y),
            'facet': FactorVector(StrVector(facet), levels=StrVector(labels))})

    cmd = 'gg <- ggplot(df, aes(x,y)) + geom_point(alpha=0.2) + ' + \
            'geom_abline(intercept=0, slope=1, color="white") + ' + \
            'stat_smooth(method="loess", color="white", size=1.5, alpha=0.2, se=FALSE) + ' + \
            'stat_smooth(method="loess", color="black") + ' + \
            'coord_cartesian(xlim=c(-0.24,6.24), ylim=c(-0.24,6.24)) + ' + \
            'xlim(0,100) + ylim(0,100) + ' + \
            'xlab("'+xlab+'") + ylab("'+ylab+'") + ' + \
            'facet_grid(facet ~ ., labeller=label_parsed)'
    R(cmd)
    if output_file:
        R('ggsave("' + output_file + '", gg, width=5, height=5)')
    else:
        print R['gg']
        raw_input('Press Enter to continue...')

def main():
    from results import get_results
    stats, args = get_results()

    plot_n_t_ratios(stats, **vars(args)) 
        
if __name__ == '__main__':
    main()

