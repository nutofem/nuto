#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"

std::vector<Eigen::VectorXd> GetTestPoints()
{
    return {Eigen::Vector3d(1. / 3., 1. / 3., 1. / 3.)};
}

BOOST_AUTO_TEST_CASE(InterpolationTetrahedronQuadratic)
{
    NuTo::Test::RunTests<NuTo::InterpolationTetrahedronQuadratic>(GetTestPoints());
}
