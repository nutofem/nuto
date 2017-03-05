#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/interpolation/Interpolation.h"


void TestN(int rDim)
{
    fakeit::Mock<NuTo::Interpolation> p;
    Method(p, GetNumNodes)       = 3;
    Method(p, GetDofDimension)   = rDim;
    Method(p, GetShapeFunctions) = Eigen::Vector3d({42, 6174, 12});

    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(rDim, 3 * rDim);
    for (int iDim = 0; iDim < rDim; ++iDim)
    {
        expected(iDim, iDim)               = 42;
        expected(iDim, iDim + rDim)        = 6174;
        expected(iDim, iDim + rDim + rDim) = 12;
    }
    BoostUnitTest::CheckEigenMatrix(p.get().GetN(Eigen::VectorXd()), expected);
}

BOOST_AUTO_TEST_CASE(InterpolationN)
{
    for (int dimension = 0; dimension < 3; ++dimension)
        TestN(dimension);
}
