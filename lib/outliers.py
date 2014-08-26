from __future__ import division

from numpy import abs, floor, array, std, mean, nan

from rpy2.robjects import DataFrame, FloatVector, StrVector, FactorVector
from rpy2.robjects.packages import importr
quantreg = importr('quantreg')

import sys

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

def chebyshev(x, y, stdevs=3.):
    stdev = std(y)
    avg = mean(y)
    thresholds = avg - stdevs * stdev, avg + stdevs * stdev
    print 'Threshold', thresholds
    out = []
    for _x, _y in zip(x, y):
        if _y < thresholds[0]  or _y > thresholds[1]:
            print 'Lost', _x, _y
            continue
        out.append((_x, _y))
    return map(array, zip(*out))

def ols_p(x, y, p_value=0.05):
    from statsmodels.api import OLS
    bonfp = OLS(x, y).fit().outlier_test()[:,2]
    for i, p in enumerate(bonfp):
        if p < p_value:
            print 'Lost', x[i], y[i]
    return zip(*[(x[i], y[i]) for i, p in enumerate(bonfp) if p >= p_value])

def main():
    pass

if __name__ == '__main__':
    sys.exit(main())
