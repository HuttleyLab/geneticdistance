from __future__ import division

from cogent import LoadTree, DNA
from cogent.evolve.models import GTR

import lib
from data import get_aln
import jsd

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from numpy import log, array

from itertools import product

import sys

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def test_shannon():
    """shannon should give shannon entropy"""
    assert_almost_equal(jsd.shannon([0,1]), 0)
    assert_almost_equal(jsd.shannon([0.25]*4), log(4))
    assert_almost_equal(jsd.shannon([1]), 0)
    assert_almost_equal(jsd.shannon([0.25,0.25,0.5]), 1.5*log(2))

def test_distribution():
    """distribution should return empirical distribution for DNA sequence"""
    al = get_aln('General', 1031).takeSeqs(('Mouse',))
    distribution = jsd.distribution(al.getSeq('Mouse'))
    st = LoadTree(tip_names=('Mouse',))
    sm = GTR()
    lf = sm.makeLikelihoodFunction(st)
    lf.setMotifProbsFromData(al)
    probs = lf.getMotifProbs()
    assert_array_almost_equal(array(probs), array(distribution))

def test_jensen_shannon():
    """jensen_shannon should calculate Jensen Shannon Divergence"""
    def jsd_alt(P, Q):
        M = (array(P) + array(Q)) / 2
        def k_l(P, Q):
            return sum(0 if p == 0. else p*log(p/q) for p, q in zip(P, Q))
        return (k_l(P, M) + k_l(Q, M))/2
   
    distributions = [([0, 1], [1, 0]),
            ([0.25]*4, [1] + [0]*3),
            ([0.2]*5, [0.2]*5),
            ([0.1, 0.3, 0.3, 0.1], [0.8, 0.1, 0.1, 0.8])]

    for d in distributions:
        assert_almost_equal(jsd.jensen_shannon(d), jsd_alt(d[0],d[1]))

    distributions += [([0, 1], [1, 0], [0.5, 0.5]),
            distributions[1] + distributions[3]]

    for d in distributions:
        assert_almost_equal(jsd.jensen_shannon(d), 
                jsd.jensen_shannon(d, map(jsd.shannon, d)))

def main():
    test_shannon()
    print 'Done test_shannon'
    test_distribution()
    print 'Done test_distribution'
    test_jensen_shannon()
    print 'Done test_jensen_shannon'
    return 0

if __name__ == '__main__':
    sys.exit(main())
