from numpy import array
import sys

from rpy2.robjects import DataFrame, FloatVector, globalenv
from rpy2.robjects import r
from rpy2.robjects.packages import importr
stats = importr('stats')

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def quantile(x, probs=None, **kw):
    return array(r.quantile(FloatVector(x), probs=FloatVector(probs), **kw))

def through_the_origin(x, y):
    df = DataFrame({'x' : FloatVector(x), 'y' : FloatVector(y)})
    s = r.summary(r.lm('y ~ 0 + x', df))
    return {'coefficient' : s.rx2('coefficients')[0],
            'stderr' : s.rx2('coefficients')[1],
            'r.squared' : s.rx2('r.squared')[0]}

def lrt_p(lhn, lha, dfn, dfa):
    LR = 2*(lha - lhn)
    df = dfa - dfn
    return 1. - r.pchisq(LR, df)[0]

def main():
    print quantile(range(10), [0.25, 0.5, 0.75])
    from numpy import arange, random
    print through_the_origin(arange(10), arange(10) + random.random(10))

if __name__ == '__main__':
    sys.exit(main())
