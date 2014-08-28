from __future__ import division

import lib

from cogent import LoadSeqs, LoadTree
from cogent import DNA

import subprocess
import csv
import tempfile

from rpy2.robjects.lib.ggplot2 import ggplot
from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r as R

from gzip import GzipFile
from cogent.evolve.pairwise_distance import ParalinearPair

from handy_r import through_the_origin

import numpy as np
import glob
import sys
import os

def get_paralinear_distances(gene, data_directory=None, third_position=False, **kw):
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
    sequences = sequences.filtered(lambda x: set(''.join(x)) <= set(DNA))

    paralinear_calc = ParalinearPair(moltype=DNA, alignment=sequences)
    paralinear_calc.run(show_progress=False)
    dists = paralinear_calc.getPairwiseDistances()

    return {frozenset(k):v for k, v in dists.items()}

def print_stats(x, y):
    stats = through_the_origin(x, y)
    print 'slope r_squared se 0.025% 0.975%'
    slope = stats['coefficient']
    se = stats['stderr']
    r_squared = stats['r.squared']
    print slope, r_squared, se, slope - 1.96*se, slope + 1.96*se
    print len(x), ' samples'

def plot_JS_EN_scatter_by_pairs(stats, output_file=None, pair=None, **kw):
    x = []
    y = []
    ya = []
    for triad in stats:
        for r in stats[triad]:
            paralinear_dists = get_paralinear_distances(r[0]['gene'], **kw)
            ns_EN = sum(r[0]['EN'][t] for t in pair)
            s_EN = sum(r[1]['EN'][t] for t in pair)
            para = paralinear_dists[pair]
            if para:
                x.append(ns_EN)
                y.append(para)
                ya.append(s_EN)
    
    print 'paralinear stats'
    print_stats(x, y)
    print 'GTR stats'
    print_stats(x, ya)
  
    df = DataFrame({'x':FloatVector(x), 'y':FloatVector(y)})
    globalenv['df'] = df
    cmd = 'gg <- ggplot(df, aes(x, y)) + geom_point(alpha=0.2) + ' + \
            'geom_abline(intercept=0, slope=1, color="white") + ' + \
            'xlab(bquote(.("'+' to '.join(pair)+'") ~ d[ENS])) + ' + \
            'ylab(bquote(.("'+' to '.join(pair)+'") ~ d[para])) + ' + \
            'coord_cartesian(xlim=c(0,1), ylim=c(0,1))'
    R(cmd)
    if output_file:
        R('ggsave("' + output_file + '", gg, width=5, height=5)')
    else:
        print R['gg']
        raw_input('Press Enter to continue...')
    return

def main():
    from results import get_results
    stats, args = get_results()

    plot_JS_EN_scatter_by_pairs(stats, **vars(args))
     
if __name__ == '__main__':
    main()

