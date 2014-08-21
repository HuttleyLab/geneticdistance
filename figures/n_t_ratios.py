from __future__ import division

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r as R

import numpy as np
import sys
import os

def plot_ENS_hexbin(stats, pair=None, output_file=None, **kw):
    pair = list(pair)
    x = [(r[0]['EN'][pair[0]]/r[0]['EN'][pair[1]])
            for t in stats for r in stats[t]]
    y = [(r[1]['EN'][pair[0]]/r[1]['EN'][pair[1]]) 
            for t in stats for r in stats[t]]

    print 'Samples', len(x)

    xlab = 'd[ENS]^{'+pair[0]+'} / d[ENS]^{'+pair[1]+'}'
    GTR = 'GTR' + ('+Gamma' if r[1]['with_rate'] else '')
    ylab = 'd['+GTR+']^{'+pair[0]+'} / d['+GTR+']^{'+pair[1]+'}'

    globalenv['df'] = DataFrame({'x': FloatVector(x), 'y': FloatVector(y)})
    cmd = 'gg <- ggplot(df, aes(x,y)) + geom_point(alpha=0.2) + ' + \
            'geom_abline(intercept=0, slope=1, color="white") + ' + \
            'stat_smooth(method="loess", color="white", size=1.5, alpha=0.2, se=FALSE) + ' + \
            'stat_smooth(method="loess", color="black") + ' + \
            'coord_cartesian(xlim=c(-0.24,6.24), ylim=c(-0.24,6.24)) + ' + \
            'xlim(0,100) + ylim(0,100) + ' + \
            'xlab(expression('+xlab+')) + ylab(expression('+ylab+'))'
    R(cmd)
    if output_file:
        R('ggsave("' + output_file + '", gg, width=5, height=5)')
    else:
        print R['gg']
        raw_input('Press Enter to continue...')

def main():
    from results import get_results
    stats, args = get_results()

    plot_ENS_hexbin(stats, **vars(args)) # vb_two
        
if __name__ == '__main__':
    main()

