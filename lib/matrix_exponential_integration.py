import cogent.maths.matrix_exponentiation as cme
from numpy.linalg import inv, eig
import numpy as np

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler', 'Von Bing Yap']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.1-dev'

class _Exponentiator(object):
    def __init__(self, Q):
        self.Q = Q

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, repr(self.Q))

class VanLoanIntegratingExponentiator(_Exponentiator):
    
    """An exponentiator that evaluates int_0^t exp(Q*s)ds * R
    using the method of Van Loan [1]. Complexity is that of the Exponentiator.
    [1] Van Loan, C. F. (1978). Computing integrals involving the matrix
    exponential. IEEE Trans. Autmat. Control 23(3), 395-404."""

    def __init__(self, Q, R=None, exponentiator=cme.RobustExponentiator):
        """
        Q -- an n x n matrix.
        R -- an n x m matrix. Defaults to the identity matrix. Can be a rank-1
        array.
        exponentiator -- Exponentiator used in Van Loan method. Defaults to
        RobustEstimator.
        """
        self.Q = Q
        Qdim = len(Q)
        if R is None:
            self.R = np.identity(Qdim)
        else:
            if len(R.shape) == 1: # Be kind to rank-1 arrays
                self.R = R.reshape((R.shape[0], 1))
            else:
                self.R = R
        Cdim = Qdim + self.R.shape[1]
        C = np.zeros((Cdim, Cdim))
        C[:Qdim,:Qdim] = Q
        C[:Qdim,Qdim:] = self.R
        self.expm = exponentiator(C)

    def __call__(self, t=1.0):
        return self.expm(t)[:len(self.Q),len(self.Q):]

class VonBingIntegratingExponentiator(_Exponentiator):

    """An exponentiator that evaluates int_0^t exp(Q*s)ds
    using the method of Von Bing."""

    def __init__(self, Q):
        """ 
        Q -- a diagonisable matrix.
        """
        self.Q = Q
        self.roots, self.evT = eig(Q)
        self.evI = inv(self.evT.T)
        # Remove following check if performance is a concern
        reQ = np.inner(self.evT*self.roots, self.evI).real
        if not np.allclose(Q, reQ): 
            raise ArithmeticError, "eigendecomposition failed"

    def __call__(self, t=1.0):
        int_roots = np.array([t if abs(x.real) < 1e-6 else
            (np.exp(x*t)-1)/x for x in self.roots])
        result = np.inner(self.evT * int_roots, self.evI)
        if result.dtype.kind == "c":
            result = np.asarray(result.real)
        result = np.maximum(result, 0.0)
        return result

