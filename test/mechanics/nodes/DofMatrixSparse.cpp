#include "BoostUnitTest.h"
#include "mechanics/nodes/DofMatrixSparse.h"
#include <sstream>

Eigen::SparseMatrix<int> Constant(int rows, int cols, int val)
{
    Eigen::SparseMatrix<int> m(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            m.coeffRef(i,j) = val;
    return m;
}

void CheckConstant(const Eigen::SparseMatrix<int>& m, int rows, int cols, int val)
{
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            BOOST_CHECK_EQUAL(m.coeff(i,j), val);
}

NuTo::DofMatrixSparse<int> Get(const NuTo::DofType& d0, const NuTo::DofType& d1)
{
    NuTo::DofMatrixSparse<int> dofMatrix;
    dofMatrix(d0, d0) = Constant(2, 2, 22);
    dofMatrix(d0, d1) = Constant(2, 3, 23);
    dofMatrix(d1, d0) = Constant(3, 2, 32);
    dofMatrix(d1, d1) = Constant(3, 3, 33);
    return dofMatrix;
}

void CheckDofMatrix(const NuTo::DofMatrixSparse<int>& m, const NuTo::DofType& d0, const NuTo::DofType& d1)
{
    CheckConstant(m(d0, d0), 2, 2, 44);
    CheckConstant(m(d0, d1), 2, 3, 46);
    CheckConstant(m(d1, d0), 3, 2, 64);
    CheckConstant(m(d1, d1), 3, 3, 66);
}

BOOST_AUTO_TEST_CASE(DofMatrixAdditionMultiplication)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);

    NuTo::DofMatrixSparse<int> dofMatrix0 = Get(dof0, dof1);
    NuTo::DofMatrixSparse<int> dofMatrix1 = Get(dof0, dof1);

    dofMatrix0 += dofMatrix1;
    CheckDofMatrix(dofMatrix0, dof0, dof1);

    dofMatrix1 *= 2;
    CheckDofMatrix(dofMatrix1, dof0, dof1);
}

BOOST_AUTO_TEST_CASE(DofMatrixUninitializedAddition)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);

    NuTo::DofMatrixSparse<int> dofMatrix0;
    NuTo::DofMatrixSparse<int> dofMatrix1 = Get(dof0, dof1);

    dofMatrix0 += dofMatrix1;
    dofMatrix0 += dofMatrix1;
    CheckDofMatrix(dofMatrix0, dof0, dof1);
}

BOOST_AUTO_TEST_CASE(DofMatrixStream)
{
    NuTo::DofType dof0("foo", 1, 0);
    NuTo::DofType dof1("bar", 1, 1);
    NuTo::DofMatrixSparse<int> dofMatrix = Get(dof0, dof1);
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

    BOOST_TEST_MESSAGE("" << dofMatrix);
}
