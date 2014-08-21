from __future__ import division

import lib

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, IntVector, StrVector, globalenv
from rpy2.robjects import r as R

import numpy as np
from scipy.stats.mstats import mquantiles

from collections import defaultdict
from results import gs_p, get_results
import sys

def get_alpha(g):
    q25, q75 = mquantiles(g, (0.25, 0.75))
    x = 2.5*q75 - 1.5*q25
    s = sum(1 for _g in g if _g < x and _g >= x-0.01)
    return 0 if s == 0 else 1/s

def plot_bar(stats, output_file=None, **kw):
    names = [r['name'] for r in stats.values()[0][0]]
    with_rates = [r['with_rate'] for r in stats.values()[0][0]]
    names = [n+('+Gamma' if w else '') for n, w in zip(names, with_rates)]

    by_dir = defaultdict(list)
    for triad in stats:
        for r in stats[triad]:
            by_dir[r[0]['from_directory']].append(r)

    for d in by_dir:
        by_dir[d] = zip(*[[gs_p(_r['gs_p']) for _r in r] for r in by_dir[d]])

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
        for j, g in enumerate(v):
            g_stats += g
            data += [dataset]*len(g)
            runs += [j]*len(g)
            print names[j], sum(1 for _g in g if _g > 0.05)/len(g)
            alpha = max(alpha, get_alpha(g))
        print 'Samples', len(g)
    labels = 'expression('+','.join(names)+')'
    
    df = DataFrame({'run':IntVector(runs), 'g_stat':FloatVector(g_stats),
        'data':StrVector(data)})
    globalenv['df'] = df
    R('library(scales)')
#            'geom_jitter(alpha=0.2, size=1) + ' + \
#            'geom_boxplot(fill=NA, outlier.size=0, size=1.5, color=alpha("white", 0.5)) + ' + \
#            'geom_boxplot(alpha=0.8, outlier.size=0) + ' + \
#            'geom_hline(yintercept=0.05, size=1.5, alpha=0.5, color="white") + ' + \
#            'geom_hline(yintercept=0.05, color="black") + ' + \
    cmd = 'gg <- ggplot(df, aes(factor(run), g_stat)) + ' + \
            'ylab("Goodness-of-Fit p-value") + xlab("Model") + ' + \
            'geom_boxplot(outlier.size=1, outlier.colour=alpha("black",'+str(alpha)+')) + ' + \
            'scale_x_discrete(labels=' + labels + ') + ' + \
            'theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + ' + \
            'facet_grid(. ~ data)'
    R(cmd)
    if output_file:
        R('ggsave("' + output_file + '", gg, width=5, height=5)')
    else:
        print R['gg']
        raw_input('Press Enter to continue...')



def main():
    stats, args = get_results()

    plot_bar(stats, **vars(args))

    return 0

if __name__ == '__main__':
    sys.exit(main())
