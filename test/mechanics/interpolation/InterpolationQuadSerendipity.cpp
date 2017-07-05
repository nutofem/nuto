#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "mechanics/interpolation/InterpolationQuadSerendipity.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    constexpr double a = 0.577350269189626;
    return {Eigen::Vector2d({-a, -a}), Eigen::Vector2d({-a, a}), Eigen::Vector2d({a, a}), Eigen::Vector2d({a, -a})};
}

BOOST_AUTO_TEST_CASE(InterpolationQuadSerendipity)
{
    NuTo::Test::RunTests<NuTo::InterpolationQuadSerendipity>(GetTestPoints());
}
