#include "BoostUnitTest.h"
#include "mechanics/dofs/DofVectorConvertEigen.h"

using namespace NuTo;

struct TestVector
{
    DofType dof0 = DofType("foo", 1);
    DofType dof1 = DofType("bar", 1);
    DofVector<double> v0;
    TestVector()
    {
        v0[dof0] = Eigen::Vector3d({1, 2, 3});
        v0[dof1] = Eigen::Vector2d({8, 9});
    }
};


BOOST_FIXTURE_TEST_CASE(DofVectorExport, TestVector)
{
    Eigen::VectorXd vExportD1 = ToEigen(v0, {dof1});
    BoostUnitTest::CheckVector(vExportD1, std::vector<double>{8, 9}, 2);

    Eigen::VectorXd vExportD0D1 = ToEigen(v0, {dof0, dof1});
    BoostUnitTest::CheckVector(vExportD0D1, std::vector<double>{1, 2, 3, 8, 9}, 5);

    Eigen::VectorXd vExportD1D0 = ToEigen(v0, {dof1, dof0});
    BoostUnitTest::CheckVector(vExportD1D0, std::vector<double>{8, 9, 1, 2, 3}, 5);
}

BOOST_FIXTURE_TEST_CASE(DofVectorImport, TestVector)
{
    Eigen::VectorXd v(5);
    v << 0, 1, 2, 3, 4;

    FromEigen(v0, v, {dof0, dof1});
    BoostUnitTest::CheckEigenMatrix(v0[dof0], Eigen::Vector3d(0, 1, 2));
    BoostUnitTest::CheckEigenMatrix(v0[dof1], Eigen::Vector2d(3, 4));
}
