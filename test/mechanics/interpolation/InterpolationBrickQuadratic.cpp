#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "nuto/mechanics/interpolation/InterpolationBrickQuadratic.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    constexpr double a = 0.577350269189626;
    return {Eigen::Vector3d({-a, -a, -a}), Eigen::Vector3d({-a, a, -a}), Eigen::Vector3d({a, a, -a}),
            Eigen::Vector3d({a, -a, -a}),  Eigen::Vector3d({-a, -a, a}), Eigen::Vector3d({-a, a, a}),
            Eigen::Vector3d({a, a, a}),    Eigen::Vector3d({a, -a, a})};
}

BOOST_AUTO_TEST_CASE(InterpolationBrickQuadratic)
{
    NuTo::Test::RunTests<NuTo::InterpolationBrickQuadratic>(GetTestPoints());
}
