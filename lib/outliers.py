from __future__ import division

from numpy import abs, floor, array, std, mean, nan

from rpy2.robjects import DataFrame, FloatVector, StrVector, FactorVector
from rpy2.robjects.packages import importr
quantreg = importr('quantreg')

import sys

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

def qlim(x, y):
    df = DataFrame({'x':FloatVector(x), 'y':FloatVector(y)})
    rq = quantreg.rq('y ~ x', df, tau=FloatVector((0.25, 0.5, 0.75)))
    print rq.rx2('coefficients')
    fv = array(rq.rx2('fitted.values'))
    return 2*max(fv[:,2] - fv[:,1]) + max(fv[:,1])

def qcrop(xlist, ylist, labels=None):
    if labels is None:
        labels = map(str, range(len(xlist)))
    x = []
    y = []
    xcrop = []
    ycrop = []
    facet = []
    for i, (onex, oney) in enumerate(zip(xlist, ylist)):
        ymax = qlim(onex, oney)
        cropx, cropy = zip(*[(nan, nan) if vy > ymax else (vx, vy)
                for vx, vy in zip(onex, oney)])
        xcrop += cropx
        ycrop += cropy
        x += onex
        y += oney
        facet += [labels[i]]*len(onex)

    df = DataFrame({'x': FloatVector(x), 'y': FloatVector(y),
            'xcrop': FloatVector(xcrop) , 'ycrop': FloatVector(ycrop),
            'facet': FactorVector(StrVector(facet), levels=StrVector(labels))})
    return df

def main():
    pass

if __name__ == '__main__':
    sys.exit(main())
