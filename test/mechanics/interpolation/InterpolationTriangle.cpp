#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "TypeTraits.h"
#include "mechanics/interpolation/InterpolationTriangle.h"


BOOST_AUTO_TEST_CASE(InterpolationTriangleCopyMove)
{
    NuTo::Test::Copy<NuTo::InterpolationTriangle>();
    NuTo::Test::Move<NuTo::InterpolationTriangle>();
}

std::vector<Eigen::VectorXd> GetTestPoints()
{
    return {Eigen::Vector2d({1. / 6., 1. / 6.}), Eigen::Vector2d({4. / 6., 1. / 6.}),
            Eigen::Vector2d({1. / 6., 4. / 6.})};
}


BOOST_AUTO_TEST_CASE(InterpolationTriangleN)
{
    NuTo::InterpolationTriangle interpolation(NuTo::eInterpolation::GAUSS, 1, 1);
    NuTo::Test::CheckShapeFunctionsAndNodePositions(interpolation);
}


BOOST_AUTO_TEST_CASE(InterpolationTriangleUnity)
{
    NuTo::InterpolationTriangle interpolation(NuTo::eInterpolation::GAUSS, 1, 1);
    NuTo::Test::CheckPartitionOfUnity(interpolation, GetTestPoints());
}
