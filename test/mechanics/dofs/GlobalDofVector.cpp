#include "BoostUnitTest.h"

#include "mechanics/dofs/GlobalDofVector.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(Access)
{
    GlobalDofVector v;

    DofType dof0("d0", 1);
    DofType dof1("d1", 1);

    v.J[dof0] = Eigen::Vector3d(0, 1, 2);
    v.K[dof0] = Eigen::Vector2d(3, 4);

    v.J[dof1] = Eigen::Vector2d(10, 11);
    v.K[dof1] = Eigen::Vector3d(12, 13, 14);

    for (int i = 0; i < 5; ++i)
    {
        BOOST_CHECK_CLOSE(v(dof0, i), i, 1.e-10);
        BOOST_CHECK_CLOSE(v(dof1, i), 10 + i, 1.e-10);
    }
}
