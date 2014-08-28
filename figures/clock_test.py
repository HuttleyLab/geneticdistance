from __future__ import division

import lib

from handy_r import lrt_p

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv, StrVector
from rpy2.robjects import r as R

from bisect import bisect
import numpy as np
import sys
import os

def plot_lrt_histograms(stats, output_file, **kw):
    gtr = []
    general = []
    gtrplusgamma = []
    for triad in stats:
        for r in stats[triad]:
            for o, l in zip((0, 2, 4), (general, gtr, gtrplusgamma)):
                p = lrt_p(r[1+o]['ll'], r[o]['ll'], r[1+o]['df'], r[o]['df'])
                l.append(p)
    
    intercepts = []
    for n, v in zip(('General', 'GTR+Gamma', 'GTR'), 
            (general, gtrplusgamma, gtr)):
        i = sum(1 for p in v if p <= 0.05)/len(v)
        intercepts.append(str(np.round(i,2)))
        print n, min(v), max(v), i, len(v)
    intercepts = 'c(' + ','.join(intercepts) + ')'

    n = len(gtr)
    globalenv['df'] = DataFrame(
            {'Model': StrVector(['general']*n + ['gtrplusgamma']*n + ['gtr']*n),
             'pvalue': FloatVector(general + gtrplusgamma + gtr)})
    cmd = 'gg <- ggplot(df, aes(pvalue, group=Model, linetype=Model)) + ' + \
        'stat_ecdf(geom="line") + ' + \
        'xlab("LRT p-value") + ylab("Empirical CDF") + ' + \
        'theme(legend.position = c(0.85, 0.15)) + ' + \
        'scale_x_continuous(breaks=c(0.05,seq(0.25,1,by=0.25)), limits=c(0,1))+'+\
        'scale_y_continuous(breaks=c(' + intercepts + \
        ', seq(0,1,by=0.25)), limits=c(0,1)) + ' + \
        'scale_linetype_discrete(labels=expression(General, GTR, GTR+Gamma))'
    R(cmd)
    if output_file:
        R('ggsave("' + output_file + '", gg, width=5, height=5)')
    else:
        print R['gg']
        raw_input('Press Enter to continue...')

def main():
    from results import get_results
    stats, args = get_results()

    plot_lrt_histograms(stats, **vars(args))
        
if __name__ == '__main__':
    main()

