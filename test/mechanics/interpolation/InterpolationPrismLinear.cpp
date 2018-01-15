#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "mechanics/interpolation/InterpolationPrismLinear.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    constexpr double a = 0.577350269189626;
    return {Eigen::Vector3d({-a, -a, -a}), Eigen::Vector3d({-a, a, -a}), Eigen::Vector3d({a, a, -a}),
            Eigen::Vector3d({a, -a, -a}),  Eigen::Vector3d({-a, -a, a}), Eigen::Vector3d({-a, a, a}),
            Eigen::Vector3d({a, a, a}),    Eigen::Vector3d({a, -a, a})};
}

BOOST_AUTO_TEST_CASE(InterpolationPrismLinear)
{
    NuTo::Test::RunTests<NuTo::InterpolationPrismLinear>(GetTestPoints());
}
