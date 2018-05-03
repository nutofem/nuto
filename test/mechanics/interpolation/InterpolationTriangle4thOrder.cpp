#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "nuto/mechanics/interpolation/InterpolationTriangle4thOrder.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    return {Eigen::Vector2d({1. / 6., 1. / 6.}), Eigen::Vector2d({4. / 6., 1. / 6.}),
            Eigen::Vector2d({1. / 6., 4. / 6.})};
}

BOOST_AUTO_TEST_CASE(InterpolationTriangle4thOrder)
{
    NuTo::InterpolationTriangle4thOrder interpolation;

    BOOST_TEST_MESSAGE("Checking copy move...");
    NuTo::Test::CheckCopyMove<NuTo::InterpolationTriangle4thOrder>();
    BOOST_TEST_MESSAGE("Checking match of shape functions and node positions...");
    NuTo::Test::CheckShapeFunctionsAndNodePositions(interpolation, 1e-8);
    BOOST_TEST_MESSAGE("Checking patrition of unity...");
    NuTo::Test::CheckPartitionOfUnity(interpolation, GetTestPoints(), 1e-8);
    BOOST_TEST_MESSAGE("Checking shape function derivatives via CDF...");
    NuTo::Test::CheckDerivativeShapeFunctionsCDF(interpolation, GetTestPoints(), 1e-9);
}
