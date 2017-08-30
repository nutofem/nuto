#include "BoostUnitTest.h"
#include "mechanics/nodes/DofMatrix.h"
#include <sstream>

NuTo::DofMatrix<int> Get(const NuTo::DofType& d0, const NuTo::DofType& d1)
{
    NuTo::DofMatrix<int> dofMatrix;
    dofMatrix(d0, d0) = Eigen::MatrixXi::Constant(2, 2, 22);
    dofMatrix(d0, d1) = Eigen::MatrixXi::Constant(2, 3, 23);
    dofMatrix(d1, d0) = Eigen::MatrixXi::Constant(3, 2, 32);
    dofMatrix(d1, d1) = Eigen::MatrixXi::Constant(3, 3, 33);
    return dofMatrix;
}

void CheckDofMatrix(const NuTo::DofMatrix<int>& m, const NuTo::DofType& d0, const NuTo::DofType& d1)
{
    BoostUnitTest::CheckEigenMatrix(m(d0, d0), Eigen::MatrixXi::Constant(2, 2, 44));
    BoostUnitTest::CheckEigenMatrix(m(d0, d1), Eigen::MatrixXi::Constant(2, 3, 46));
    BoostUnitTest::CheckEigenMatrix(m(d1, d0), Eigen::MatrixXi::Constant(3, 2, 64));
    BoostUnitTest::CheckEigenMatrix(m(d1, d1), Eigen::MatrixXi::Constant(3, 3, 66));
}

BOOST_AUTO_TEST_CASE(DofMatrixAdditionMultiplication)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);

    NuTo::DofMatrix<int> dofMatrix0 = Get(dof0, dof1);
    NuTo::DofMatrix<int> dofMatrix1 = Get(dof0, dof1);

    dofMatrix0 += dofMatrix1;
    CheckDofMatrix(dofMatrix0, dof0, dof1);

    dofMatrix1 *= 2;
    CheckDofMatrix(dofMatrix1, dof0, dof1);
}

BOOST_AUTO_TEST_CASE(DofMatrixUninitializedAddition)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);

    NuTo::DofMatrix<int> dofMatrix0;
    NuTo::DofMatrix<int> dofMatrix1 = Get(dof0, dof1);

    dofMatrix0 += dofMatrix1;
    dofMatrix0 += dofMatrix1;
    CheckDofMatrix(dofMatrix0, dof0, dof1);
}

BOOST_AUTO_TEST_CASE(DofVectorStream)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);
    NuTo::DofMatrix<int> dofMatrix = Get(dof0, dof1);
    std::stringstream ss;
    ss << dofMatrix;
    // header:
    BOOST_CHECK(ss.str().find("=== 0 0 ===") != std::string::npos);
    BOOST_CHECK(ss.str().find("=== 0 1 ===") != std::string::npos);
    BOOST_CHECK(ss.str().find("=== 1 0 ===") != std::string::npos);
    BOOST_CHECK(ss.str().find("=== 1 1 ===") != std::string::npos);
    // values:
    BOOST_CHECK(ss.str().find("22") != std::string::npos);
    BOOST_CHECK(ss.str().find("23") != std::string::npos);
    BOOST_CHECK(ss.str().find("32") != std::string::npos);
    BOOST_CHECK(ss.str().find("33") != std::string::npos);
}
