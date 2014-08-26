from __future__ import division

import lib

from outliers import qcrop

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
quantreg = importr('quantreg')

import numpy as np
import sys

def label(row):
    if row['with_rate']:
        return row['name']+'+Gamma'
    return row['name']

def plot_JS_EN_scatter_by_pairs(stats, output_file=None, outliers=0., **kw):
    xlist = []
    ylist = []
    labellist = []
    for o in range(0,len(stats.values()[0][0]),2):
        labellist.append(label(stats.values()[0][0][o+1]))
        x = []
        y = []
        for triad in stats:
            for r in stats[triad]:
                    max_jsd = -1, None
                    for pair in r[o+1]['js']:
                        if pair == triad:
                            continue
                        jsd = r[o+1]['js'][pair]
                        if jsd > max_jsd[0]:
                            ns_EN = sum(r[o]['EN'][t] for t in pair)
                            s_EN = sum(r[o+1]['EN'][t] for t in pair)
                            max_jsd = jsd, s_EN - ns_EN
                    x.append(max_jsd[0])
                    y.append(max_jsd[1])
        xlist.append(x)
        ylist.append(y)

    print 'Samples', len(xlist[0])

    xlabel = 'Max Pairwise JSD Over Triad'
    ylabel = 'Corresponding Genetic Distance Error'

    globalenv['df'] = qcrop(xlist, ylist, labellist)

    cmd = 'gg <- ggplot(df, aes(x,y)) + ' + \
            'geom_point(aes(xcrop, ycrop), alpha=0.2) + ' + \
            'xlab("' + xlabel + '") + ylab("' + ylabel + '") + ' + \
            'stat_quantile(colour="white", size=1.5, alpha=0.5) +  ' + \
            'stat_quantile(colour="black") +  ' + \
            'facet_grid(facet ~ ., space="free", scale="free", labeller=label_parsed)'
    R(cmd)
    if output_file:
        R('ggsave("' + output_file + '", gg, width=5, height=5)')
    else:
        print R['gg']
        raw_input('Press Enter to continue...')
    
def main():
    from results import get_results
    stats, args = get_results()

    plot_JS_EN_scatter_by_pairs(stats, **vars(args))

    return 0
     
if __name__ == '__main__':
    main()

