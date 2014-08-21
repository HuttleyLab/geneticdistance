import lib
import matrix_exponential_integration as expm
from numpy.testing import assert_almost_equal, assert_raises, \
    assert_array_almost_equal
from numpy import array, dot, diag, exp
import cogent.maths.matrix_exponentiation as cmme

def test_van_loan_integrating_exponentiator():
    """VanLoanIntegratingExponentiator should reproduce Felsenstein
    analytic result, should throw if we pass it a defected matrix and ask
    it to use CheckedExponentiator, will work with a defective matrix (that
    we can integrate by hand) if we use the default RobustExponentiator,
    and should work for different choices of R and exponentiatior."""
    # Result from Von Bing's R code
    result = 0.7295333
    q = array([[0.5, 0.2, 0.1, 0.2]]*4)
    for i in range(4):
        q[i, i] = 0.
        q[i, i] = -sum(q[i,])
    p0 = array([0.2, 0.3, 0.3, 0.2])

    I = expm.VanLoanIntegratingExponentiator(q, -diag(q))(1.0)
    assert_almost_equal(dot(p0, I), result)

    assert_raises(ArithmeticError,
            expm.VanLoanIntegratingExponentiator, 
            q, -diag(q), cmme.CheckedExponentiator)

    Q = array([[1., 1.], [0., 1.]])
    def integral(t):
        return array([[exp(t)-1., exp(t)*(t-1.)+1.], [0., exp(t)-1.]])
    assert_almost_equal(expm.VanLoanIntegratingExponentiator(Q)(1.),
        integral(1.))
    assert_almost_equal(expm.VanLoanIntegratingExponentiator(Q)(2.),
        integral(2.))

    R = array([[1.],[1.]])
    assert_almost_equal(expm.VanLoanIntegratingExponentiator(Q, R,
        cmme.TaylorExponentiator)(1.), dot(integral(1.), R))

def test_von_bing_integrating_exponentiator():
    """VonBingIntegratingExponentiator should reproduce Felsenstein 
    analytic result, should throw if we pass it a defective matrix, and
    should match results obtained from VanLoanIntegratingExponentiator for
    a diagonisable matrix."""
    # Result from Von Bing's R code.
    result = 0.7295333
    q = array([[0.5, 0.2, 0.1, 0.2]]*4)
    for i in range(4):
        q[i, i] = 0.
        q[i, i] = -sum(q[i,])
    p0 = array([0.2, 0.3, 0.3, 0.2])

    I = expm.VonBingIntegratingExponentiator(q)(1.0)
    assert_almost_equal(dot(dot(p0, I), -diag(q)), result)

    assert_raises(ArithmeticError,
            expm.VonBingIntegratingExponentiator,
            array([[1., 1.], [0., 1.]]))

    p = array([[ 0.86758487,  0.05575623,  0.0196798 ,  0.0569791 ],
    [ 0.01827347,  0.93312148,  0.02109664,  0.02750842],
    [ 0.04782582,  0.1375742 ,  0.80046869,  0.01413129],
    [ 0.23022035,  0.22306947,  0.06995306,  0.47675713]])

    assert_array_almost_equal(expm.VonBingIntegratingExponentiator(p)(1.),
            expm.VanLoanIntegratingExponentiator(p,
                exponentiator=cmme.FastExponentiator)(1.))
    assert_array_almost_equal(expm.VonBingIntegratingExponentiator(p)(2.),
            expm.VanLoanIntegratingExponentiator(p, 
                exponentiator=cmme.FastExponentiator)(2.))



