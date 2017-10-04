#include "BoostUnitTest.h"
#include "mechanics/dofs/DofVector.h"
#include <sstream>

BOOST_AUTO_TEST_CASE(DofVectorAddition)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);

    NuTo::DofVector<double> dofVector0;
    NuTo::DofVector<double> dofVector1;

    dofVector0[dof0] = Eigen::Vector3d({1, 2, 3});
    dofVector1[dof0] = Eigen::Vector3d({10, 20, 30});

    dofVector0[dof1] = Eigen::Vector2d({8, 9});
    dofVector1[dof1] = Eigen::Vector2d({80, 90});

    dofVector0 += dofVector1;

    BoostUnitTest::CheckVector(dofVector0[dof0], std::vector<double>({11, 22, 33}), 3);
    BoostUnitTest::CheckVector(dofVector0[dof1], std::vector<double>({88, 99}), 2);
}

BOOST_AUTO_TEST_CASE(DofVectorScalarMultiplication)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);

    NuTo::DofVector<double> dofVector0;
    dofVector0[dof0] = Eigen::Vector3d({1, 2, 3});
    dofVector0[dof1] = Eigen::Vector2d({8, 9});

    dofVector0 *= 2.;
    BoostUnitTest::CheckVector(dofVector0[dof0], std::vector<double>({2, 4, 6}), 3);
    BoostUnitTest::CheckVector(dofVector0[dof1], std::vector<double>({16, 18}), 2);

    NuTo::DofVector<double> dofVector1 = dofVector0 * 0.5;
    BoostUnitTest::CheckVector(dofVector1[dof0], std::vector<double>({1, 2, 3}), 3);
    BoostUnitTest::CheckVector(dofVector1[dof1], std::vector<double>({8, 9}), 2);
}

BOOST_AUTO_TEST_CASE(DofVectorUninitializedAddition)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);

    NuTo::DofVector<double> dofVector0;
    NuTo::DofVector<double> dofVector1;

    dofVector1[dof0] = Eigen::Vector3d({11, 22, 33});
    dofVector1[dof1] = Eigen::Vector2d({88, 99});

    dofVector0 += dofVector1;

    BoostUnitTest::CheckVector(dofVector0[dof0], std::vector<double>({11, 22, 33}), 3);
    BoostUnitTest::CheckVector(dofVector0[dof1], std::vector<double>({88, 99}), 2);
}

BOOST_AUTO_TEST_CASE(DofVectorStream)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);
    NuTo::DofVector<double> dofVector;
    dofVector[dof0] = Eigen::Vector3d::Zero();
    dofVector[dof1] = Eigen::Vector3d::Zero();
    std::stringstream ss;
    ss << dofVector;
}

