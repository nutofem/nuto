#include "BoostUnitTest.h"

#include "mechanics/dofs/GlobalDofVector.h"

using namespace NuTo;

struct TestVector
{
    GlobalDofVector v;

    DofType dof0 = DofType("d0", 1);
    DofType dof1 = DofType("d1", 1);

    DofVector<double> v0;
    DofVector<double> v1;

    TestVector()
    {
        v.J[dof0] = Eigen::Vector3d(0, 1, 2);
        v.K[dof0] = Eigen::Vector2d(3, 4);

        v.J[dof1] = Eigen::Vector2d(10, 11);
        v.K[dof1] = Eigen::Vector3d(12, 13, 14);
    }
};

BOOST_FIXTURE_TEST_CASE(Access, TestVector)
{
    for (int i = 0; i < 5; ++i)
    {
        BOOST_CHECK_CLOSE(v(dof0, i), i, 1.e-10);
        BOOST_CHECK_CLOSE(v(dof1, i), 10 + i, 1.e-10);
    }
}

BOOST_FIXTURE_TEST_CASE(Calc, TestVector)
{
    v += v;
    for (int i = 0; i < 5; ++i)
    {
        BOOST_CHECK_CLOSE(v(dof0, i), 2 * i, 1.e-10);
        BOOST_CHECK_CLOSE(v(dof1, i), 2 * (10 + i), 1.e-10);
    }

    v -= v;

    for (int i = 0; i < 5; ++i)
    {
        BOOST_CHECK_SMALL(v(dof0, i), 1.e-10);
        BOOST_CHECK_SMALL(v(dof1, i), 1.e-10);
    }
}
