#include "BoostUnitTest.h"
#include "mechanics/nodes/DofVector.h"
#include <eigen3/Eigen/Core>

BOOST_AUTO_TEST_CASE(DofVectorAddition)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);

    NuTo::DofVector<Eigen::VectorXd> dofVector0;
    NuTo::DofVector<Eigen::VectorXd> dofVector1;

    dofVector0[dof0] = Eigen::Vector3d({1, 2, 3});
    dofVector1[dof0] = Eigen::Vector3d({10, 20, 30});

    dofVector0[dof1] = Eigen::Vector2d({8, 9});
    dofVector1[dof1] = Eigen::Vector2d({80, 90});

    dofVector0 += dofVector1;

    BoostUnitTest::CheckVector(dofVector0[dof0], std::vector<double>({11, 22, 33}), 3);
    BoostUnitTest::CheckVector(dofVector0[dof1], std::vector<double>({88, 99}), 2);
}

BOOST_AUTO_TEST_CASE(DofVectorUninitializedAddition)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);

    NuTo::DofVector<Eigen::VectorXd> dofVector0;
    NuTo::DofVector<Eigen::VectorXd> dofVector1;

    dofVector1[dof0] = Eigen::Vector3d({11, 22, 33});
    dofVector1[dof1] = Eigen::Vector2d({88, 99});

    dofVector0 += dofVector1;

    BoostUnitTest::CheckVector(dofVector0[dof0], std::vector<double>({11, 22, 33}), 3);
    BoostUnitTest::CheckVector(dofVector0[dof1], std::vector<double>({88, 99}), 2);
}
