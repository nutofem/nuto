#include "BoostUnitTest.h"
#include "nuto/mechanics/cell/Matrix.h"


void TestN(int dim)
{
    int numNodes = 3;
    Eigen::Vector3d shapeFunctions({42, 6174, 12});

    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(dim, 3 * dim);
    for (int iDim = 0; iDim < dim; ++iDim)
    {
        expected(iDim, iDim) = 42;
        expected(iDim, iDim + dim) = 6174;
        expected(iDim, iDim + dim + dim) = 12;
    }
    BoostUnitTest::CheckEigenMatrix(NuTo::Matrix::N(shapeFunctions, numNodes, dim), expected);
}

BOOST_AUTO_TEST_CASE(InterpolationN)
{
    TestN(1);
    TestN(2);
    TestN(3);
}
