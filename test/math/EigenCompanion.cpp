#include "BoostUnitTest.h"
#include "math/EigenCompanion.h"

using namespace NuTo::EigenCompanion;

BOOST_AUTO_TEST_CASE(InitList)
{
    BoostUnitTest::CheckVector(ToEigen({5, 4, 3, 2, 1}), std::vector<double>{5, 4, 3, 2, 1}, 5);
}

BOOST_AUTO_TEST_CASE(To42D)
{
    BoostUnitTest::CheckEigenMatrix(To3D(ToEigen(42)), Eigen::Vector3d(42, 0, 0));
    BoostUnitTest::CheckEigenMatrix(To3D(ToEigen({1, 2})), Eigen::Vector3d(1, 2, 0));
    BoostUnitTest::CheckEigenMatrix(To3D(ToEigen({1, 2, 3})), Eigen::Vector3d(1, 2, 3));
}
