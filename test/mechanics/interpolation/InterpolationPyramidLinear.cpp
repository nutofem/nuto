#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "mechanics/interpolation/InterpolationPyramidLinear.h"

// Test points should be inside the pyramid element
// The base is a quad in xy plane (z=0) [-1,1] * [-1,1]
// and the top is at (0,0,1)
std::vector<Eigen::VectorXd> GetTestPoints()
{
    return {Eigen::Vector3d({-0.7, -0.7, 0}), Eigen::Vector3d({0.7, -0.7, 0}),    Eigen::Vector3d({0.7, 0.7, 0}),
            Eigen::Vector3d({-0.7, 0.7, 0}),  Eigen::Vector3d({-0.1, -0.1, 0.8}), Eigen::Vector3d({0.1, -0.1, 0.8}),
            Eigen::Vector3d({0.1, 0.1, 0.8}), Eigen::Vector3d({-0.1, 0.1, 0.8}),  Eigen::Vector3d({0.0, 0.0, 0.99})};
}

BOOST_AUTO_TEST_CASE(InterpolationPyramidLinear)
{
    NuTo::Test::RunTests<NuTo::InterpolationPyramidLinear>(GetTestPoints());
}
