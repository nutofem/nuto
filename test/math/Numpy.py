import nuto
import numpy as np
import sys
from ctypes import c_int



def test_double():
    denseTest = nuto.NumpyDense
    toCpp = np.array([[1., 2.],[3., 4.],[5., 6.]])
    expected = np.array([[2., 4.],[6., 8.],[10., 12.]])

    fromCpp = denseTest.CheckMatrixXd(toCpp)
    if np.max(np.abs(fromCpp - expected)) > 1.e-10:
        print "expected \n", expected
        print "got \n", fromCpp
        sys.exit(-1)

def test_int():
    denseTest = nuto.NumpyDense
    toCpp = np.array([[1, 2],[3, 4],[5, 6]], dtype=c_int)
    expected = np.array([[2, 4],[6, 8],[10, 12]])

    fromCpp = denseTest.CheckMatrixXi(toCpp)
    if np.max(np.abs(fromCpp - expected)) > 1.e-10:
        print "expected \n", expected
        print "got \n", fromCpp
        sys.exit(-1)

test_double()
test_int()