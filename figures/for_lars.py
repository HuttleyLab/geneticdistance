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
    total_sig = 0
    sig_with_human_gt_mouse = 0
    samples = 0
    for triad in stats:
        for r in stats[triad]:
            p = lrt_p(r[1]['ll'], r[0]['ll'], r[1]['df'], r[0]['df'])
            if p <= 0.05:
                total_sig += 1
                if r[0]['EN']['Human'] > r[0]['EN']['Mouse']:
                    sig_with_human_gt_mouse += 1
            samples += 1
    print samples, 'samples'
    print total_sig, 'significant'
    print sig_with_human_gt_mouse, 'ENS_Human > ENS_Mouse'

def main():
    from results import get_results
    stats, args = get_results()

    plot_lrt_histograms(stats, **vars(args))
        
if __name__ == '__main__':
    main()

