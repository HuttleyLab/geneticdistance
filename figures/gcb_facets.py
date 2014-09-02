from __future__ import division

import lib

from outliers import qcrop
from handy_r import quantile

from cogent import LoadSeqs, LoadTree
from cogent import DNA

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r as R

from gzip import GzipFile

import itertools
import numpy as np
import glob
import sys
import os

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

def get_sequences(gene, data_directory=None, third_position=False, **kw):
    filenames = glob.glob(os.path.join(data_directory, gene+'.fasta*'))
    assert len(filenames) == 1, 'Wrong number of alignment files for ' + gene
    filename = filenames[0]
    if filename.endswith('.fasta'):
        with open(filename) as fastafile:
            fastadata = fastafile.read()
    elif filename.endswith('.fasta.gz'):
        with GzipFile(filename) as fastafile:
            fastadata = fastafile.read()
    else:
        raise RuntimeError(gene + ' file could not be read')

    sequences = LoadSeqs(data=fastadata)
    if third_position:
        indices = [(i, i+1) for i in range(len(sequences))[2::3]]
        pos3 = sequences.addFeature('pos3', 'pos3', indices)
        sequences = pos3.getSlice()

    return sequences

def get_gc_biases(sequences, triad):
    gcs = {}
    sequences = sequences.takeSeqs(triad)
    sequences = sequences.filtered(lambda x: set(''.join(x)) <= set(DNA))
    for name in triad:
        seq = sequences.getGappedSeq(name)
        num_gc = sum(1 for p in seq if p in 'GC')
        gcs[name] = num_gc/len(seq)

    pairwise_gc_diffs = {}
    for n1, n2 in itertools.combinations(triad, 2):
        pairwise_gc_diffs[frozenset([n1, n2])] = abs(gcs[n1] - gcs[n2])

    return pairwise_gc_diffs

def plot_JS_EN_scatter_by_pairs(stats, output_file=None, **kw):
    xlist = []
    ylist = []
    labellist = []
    for o in range(0,len(stats.values()[0][0]),2):
        labellist.append(label(stats.values()[0][0][o+1]))
        x = []
        y = []
        gene = None
        for i, triad in enumerate(stats):
            for j, r in enumerate(stats[triad]):
                if gene != r[o]['gene']:
                    gene = r[o]['gene'].split('_')[0]
                    sequences = get_sequences(gene, **kw)
                gc_diffs = get_gc_biases(sequences, triad)
                pair = sorted([(v, k) for k, v in gc_diffs.items()]).pop()[1]
                ns_EN_alt = sum(r[o]['EN'][t] for t in pair)
                ns_EN_null = sum(r[o+1]['EN'][t] for t in pair)
                y.append(ns_EN_null / ns_EN_alt)
                x.append(gc_diffs[pair])
        xlist.append(x)
        ylist.append(y)
    
        print labellist[-1], str(len(x)), 'samples' 
        n20 = [_y for _x, _y in zip(x,y) if _x >= 0.19 and _x <= 0.21]
        if len(n20) > 0:
            print 'mean', 'quantile025', 'quantile975', 'length'
            print np.mean(n20), quantile(n20, [0.025,0.975]), len(n20)
    
    globalenv['df'] = qcrop(xlist, ylist, labellist)

    xlabel = 'Max Pairwise Change in G+C Bias Over Triad'
    ylabel = 'Corresponding Proportional Genetic Distance Error'
    cmd = 'gg <- ggplot(df, aes(x,y)) + ' + \
            'geom_point(aes(xcrop, ycrop), alpha=0.2) + ' + \
            'stat_smooth(method="loess", color="white", size=1.5, alpha=0.2, se=FALSE) + ' + \
            'stat_smooth(method="loess", color="black") + ' + \
            'xlab("'+xlabel+'") + ylab("'+ylabel+'") + ' + \
            'facet_grid(facet ~ ., space="free", scale="free", labeller=label_parsed)'
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

