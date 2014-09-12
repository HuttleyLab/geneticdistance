import numpy as np
from scipy_optimize import newton

from cogent.evolve.substitution_model import General
from cogent.evolve.substitution_calculation import (
        CalcDefn, NonParamDefn, ExpDefn)
from cogent.maths.matrix_exponentiation import EigenExponentiator
from cogent.maths.matrix_exponential_integration import (
        VanLoanIntegratingExponentiator, VonBingIntegratingExponentiator)

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.3-dev'

class CalcQd(object):
    def __init__(self, exp, calcExMat, word_probs, mprobs_matrix, *params):
        self.mprobs_matrix = mprobs_matrix
        self.R = calcExMat(word_probs, *params)
        self.R *= mprobs_matrix
        row_totals = self.R.sum(axis=1)
        self.R -= np.diag(row_totals)
        self.guess_alpha = 1. / (word_probs * row_totals).sum()
        self.Rd = exp(self.R)
        if isinstance(self.Rd, EigenExponentiator):
            iexp = VonBingIntegratingExponentiator
            iexpR = iexp(self.R)
            self.iRd = lambda t: np.dot(iexpR(t), -np.diag(self.R))
        else:
            iexp = VanLoanIntegratingExponentiator
            self.iRd = iexp(self.R, -np.diag(self.R), exp)
        self._fprimealpha = None

    def _f(self, alpha):
        ENS = np.dot(self.mprobs_matrix, self.iRd(alpha*self.distance))
        f = ENS - self.distance
        return f
    
    def _fprime(self, alpha):
        df = np.dot(self.mprobs_matrix, self._fprimecommon(alpha))
        return df
    
    def _fprime2(self, alpha):
        d2f = self.distance*np.dot(np.dot(self.mprobs_matrix, self.R),
                self._fprimecommon(alpha))
        return d2f

    def _fprimecommon(self, alpha):
        if self._fprimealpha != alpha:
            self._fprimevalue = -self.distance * \
                    np.dot(self.Rd(alpha * self.distance), np.diag(self.R))
            self._fprimealpha = alpha
        return self._fprimevalue

    def getQ(self, distance):
        self.distance = distance
        try:
            self.alpha = newton(self._f, self.guess_alpha, self._fprime,
                    fprime2=self._fprime2)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.alpha = self.guess_alpha
        Q = self.alpha * self.R
        return Q

    def getP(self):
        return self.Rd(self.alpha * self.distance)


class GeneralBen(General):
    """A continuous substitution model with one free parameter for each and 
    every possible instantaneous substitution for which branch length equals
    expected number of substitutions along that branch"""
    
    def __init__(self, *args, **kwargs):
        if 'name' not in kwargs:
            kwargs['name'] = 'GeneralBen'
        super(GeneralBen, self).__init__(*args, **kwargs)

    def calcQ(self):
        raise AttributeError("'GeneralBen' object has no attribute 'calcQ'")

    def makeQdDefn(self, word_probs, mprobs_matrix, rate_params):
        expm = NonParamDefn('expm')
        exp = ExpDefn(expm)
        calcExMat = CalcDefn(lambda: self.calcExchangeabilityMatrix,
            name='calculate_exchangeability_matrix')()
        Qd = CalcDefn(CalcQd, name='Qd')(exp, calcExMat, word_probs, 
                mprobs_matrix, *rate_params)
        return Qd

    def makeContinuousPsubDefn(self, word_probs, mprobs_matrix, distance, rate_params):
        Qd = self.makeQdDefn(word_probs, mprobs_matrix, rate_params)
        Q = CalcDefn(lambda cQd, t: cQd.getQ(t), name='Q')(Qd, distance)
        P = CalcDefn(lambda cQd, Q: cQd.getP(), name='psubs')(Qd, Q)
        return P




