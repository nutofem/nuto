#include "BoostUnitTest.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include <sstream>

using namespace NuTo;

struct TestVectors
{
    DofType dof0 = DofType("foo", 1);
    DofType dof1 = DofType("bar", 1);

    DofVector<double> v0;
    DofVector<double> v1;

    TestVectors()
    {
        v0[dof0] = Eigen::Vector3d({1, 2, 3});
        v0[dof1] = Eigen::Vector2d({8, 9});

        v1[dof0] = Eigen::Vector3d({10, 20, 30});
        v1[dof1] = Eigen::Vector2d({80, 90});
    }
};

BOOST_FIXTURE_TEST_CASE(DofVectorAddition, TestVectors)
{
    v0 += v1;
    BoostUnitTest::CheckEigenMatrix(v0[dof0], Eigen::Vector3d(11, 22, 33));
    BoostUnitTest::CheckEigenMatrix(v0[dof1], Eigen::Vector2d(88, 99));
    BoostUnitTest::CheckEigenMatrix(v1[dof0], Eigen::Vector3d(10, 20, 30));
    BoostUnitTest::CheckEigenMatrix(v1[dof1], Eigen::Vector2d(80, 90));
}

BOOST_FIXTURE_TEST_CASE(DofVectorScalarMultiplication, TestVectors)
{
    v0 *= 2.;
    BoostUnitTest::CheckEigenMatrix(v0[dof0], Eigen::Vector3d(2, 4, 6));
    BoostUnitTest::CheckEigenMatrix(v0[dof1], Eigen::Vector2d(16, 18));

    DofVector<double> v = v0 * 0.5;
    BoostUnitTest::CheckEigenMatrix(v[dof0], Eigen::Vector3d(1, 2, 3));
    BoostUnitTest::CheckEigenMatrix(v[dof1], Eigen::Vector2d(8, 9));
}

BOOST_FIXTURE_TEST_CASE(DofVectorUninitializedAddition, TestVectors)
{
    DofVector<double> v;
    v += v0 + v1;
    BoostUnitTest::CheckEigenMatrix(v[dof0], Eigen::Vector3d(11, 22, 33));
    BoostUnitTest::CheckEigenMatrix(v[dof1], Eigen::Vector2d(88, 99));
}

BOOST_FIXTURE_TEST_CASE(DofVectorStream, TestVectors)
{
    std::stringstream ss;
    ss << v0;
}

BOOST_FIXTURE_TEST_CASE(DofTypeAccess, TestVectors)
{
    auto dofTypes = v0.DofTypes();
    BOOST_CHECK_EQUAL(dofTypes.size(), 2);
    std::sort(dofTypes.begin(), dofTypes.end(), CompareDofType());
    BOOST_CHECK_EQUAL(dofTypes[0].Id(), dof0.Id());
    BOOST_CHECK_EQUAL(dofTypes[1].Id(), dof1.Id());
}
