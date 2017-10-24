#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    return {Eigen::VectorXd::Zero(1)};
}

BOOST_AUTO_TEST_CASE(InterpolationTrussLinear)
{
    NuTo::Test::RunTests<NuTo::InterpolationTrussLinear>(GetTestPoints());
}
