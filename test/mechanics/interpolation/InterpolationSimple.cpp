#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/interpolation/InterpolationSimple.h"


void TestN(int rDim)
{
    fakeit::Mock<NuTo::InterpolationSimple> p;
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
    BoostUnitTest::CheckEigenMatrix(p.get().GetN(Eigen::Vector2d({0.0, 0.0})), expected);
}

BOOST_AUTO_TEST_CASE(InterpolationN)
{
    for (int dimension = 1; dimension < 3; ++dimension)
        TestN(dimension);
}
