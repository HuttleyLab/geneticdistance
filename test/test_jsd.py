from __future__ import division

from cogent import LoadTree, DNA
from cogent.evolve.models import GTR

import lib
import jsd

from numpy.testing import assert_almost_equal, assert_array_almost_equal
from numpy import log, array

from itertools import product

def test_shannon():
    """shannon should give shannon entropy"""
    assert_almost_equal(jsd.shannon([0,1]), 0)
    assert_almost_equal(jsd.shannon([0.25]*4), log(4))
    assert_almost_equal(jsd.shannon([1]), 0)
    assert_almost_equal(jsd.shannon([0.25,0.25,0.5]), 1.5*log(2))

def test_distribution():
    """distribution should return empirical distribution for DNA sequence"""
    st = LoadTree(tip_names='a')
    sm = GTR()
    lf = sm.makeLikelihoodFunction(st)
    al = lf.simulateAlignment(1000)
    lf.setMotifProbsFromData(al)
    probs = lf.getMotifProbs()
    distribution = jsd.distribution(al.getSeq('a'))
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

def test_column_dist():
    """column_dist should calculate empirical column distribution"""
    st = LoadTree(tip_names='abc')
    sm = GTR()
    lf = sm.makeLikelihoodFunction(st)
    al = lf.simulateAlignment(1000)
    freqs = dict(zip([''.join(w) for w in product(*[DNA.Alphabet]*3)], [0]*64))
    for c in [''.join(w) for w in zip(*al.todict().values())]:
        freqs[c] += 1
    for k in freqs:
        freqs[k] /= 1000
    dist = jsd.column_dist(al)
    for k in set(freqs.keys()).union(dist.keys()):
        assert_almost_equal(dist[k], freqs[k])
