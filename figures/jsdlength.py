from __future__ import division

import sys
import os

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
quantreg = importr('quantreg')

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

def plot_JS_EN_scatter_by_pairs(stats, output_file=None, pair=None, **kw):
    x = []
    y = []
    for triad in stats:
        for r in stats[triad]:
    #        pair = sorted([(v, k) 
    #            for k, v in r[0]['js'].items() if len(k)==2]).pop()[1]
            x.append(r[0]['js'][pair])
            y.append(sum(r[0]['EN'][t] for t in pair))

    title = str(len(x)) + ' samples'
    if output_file:
        title = output_file + ', ' + title
    print title

    globalenv['df'] = qcrop([x], [y])

    cmd = 'gg <- ggplot(df, aes(x,y)) + ' + \
            'geom_point(aes(xcrop, ycrop), alpha=0.2) + ' + \
            'stat_smooth(method="loess", color="white", size=1.5, alpha=0.2, se=FALSE) + ' + \
            'stat_smooth(method="loess", color="black") + ' + \
            'xlab("'+' to '.join(pair)+' JSD") + ' + \
            'ylab(bquote(.("'+' to '.join(pair)+'") ~ d[ENS])) + coord_flip()'
    R(cmd)
    if output_file:
        R('ggsave("'+output_file+'", gg, width=5, height=5)')
    else:
        print R['gg']
        raw_input('Press Enter to continue...')

def main():
    from results import get_results
    stats, args = get_results()

    plot_JS_EN_scatter_by_pairs(stats, **vars(args))
     
if __name__ == '__main__':
    main()

