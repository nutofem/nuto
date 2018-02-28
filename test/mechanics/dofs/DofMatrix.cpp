#include "BoostUnitTest.h"
#include "mechanics/dofs/DofMatrix.h"
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
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);

    NuTo::DofMatrix<int> dofMatrix0 = Get(dof0, dof1);
    NuTo::DofMatrix<int> dofMatrix1 = Get(dof0, dof1);

    dofMatrix0 += dofMatrix1;
    CheckDofMatrix(dofMatrix0, dof0, dof1);

    dofMatrix1 *= 2;
    CheckDofMatrix(dofMatrix1, dof0, dof1);
}

BOOST_AUTO_TEST_CASE(DofMatrixUninitializedAddition)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);

    NuTo::DofMatrix<int> dofMatrix0;
    NuTo::DofMatrix<int> dofMatrix1 = Get(dof0, dof1);

    dofMatrix0 += dofMatrix1;
    dofMatrix0 += dofMatrix1;
    CheckDofMatrix(dofMatrix0, dof0, dof1);
}

BOOST_AUTO_TEST_CASE(DofMatrixStream)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);
    NuTo::DofMatrix<int> dofMatrix = Get(dof0, dof1);
    std::stringstream ss;
    ss << dofMatrix;

    BOOST_CHECK(ss.str().find("22") != std::string::npos);
    BOOST_CHECK(ss.str().find("23") != std::string::npos);
    BOOST_CHECK(ss.str().find("32") != std::string::npos);
    BOOST_CHECK(ss.str().find("33") != std::string::npos);

    BOOST_TEST_MESSAGE("" << dofMatrix);
}

BOOST_AUTO_TEST_CASE(DofTypeAccess)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);
    NuTo::DofMatrix<int> dofMatrix = Get(dof0, dof1);
    auto dofTypes = dofMatrix.DofTypes();

    BOOST_CHECK_EQUAL(dofTypes.size(), 2);
    std::sort(dofTypes.begin(), dofTypes.end(), NuTo::CompareDofType());
    BOOST_CHECK_EQUAL(dofTypes[0].Id(), dof0.Id());
    BOOST_CHECK_EQUAL(dofTypes[1].Id(), dof1.Id());
}

BOOST_AUTO_TEST_CASE(UnusualDofTypeAccess)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);
    NuTo::DofMatrix<int> partiallyFilledDofMatrix;
    partiallyFilledDofMatrix(dof1, dof0) = Eigen::Matrix2i::Zero();
    auto dofTypes = partiallyFilledDofMatrix.DofTypes();

    BOOST_CHECK_EQUAL(dofTypes.size(), 2);
    std::sort(dofTypes.begin(), dofTypes.end(), NuTo::CompareDofType());
    BOOST_CHECK_EQUAL(dofTypes[0].Id(), dof0.Id());
    BOOST_CHECK_EQUAL(dofTypes[1].Id(), dof1.Id());
}
