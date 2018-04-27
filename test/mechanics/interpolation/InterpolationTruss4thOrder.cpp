#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "nuto/mechanics/interpolation/InterpolationTruss4thOrder.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    return {Eigen::VectorXd::Zero(1), Eigen::VectorXd::Constant(1, -0.89), Eigen::VectorXd::Constant(1, 0.89)};
}

BOOST_AUTO_TEST_CASE(InterpolationTruss4thOrder)
{
    NuTo::Test::RunTests<NuTo::InterpolationTruss4thOrder>(GetTestPoints());
}
