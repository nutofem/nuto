#include "BoostUnitTest.h"
#include "mechanics/dofs/DofMatrixSparseConvertEigen.h"

BOOST_AUTO_TEST_CASE(DofMatrixExportEigen)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);

    NuTo::DofMatrixSparse<double> m;
    m(dof0, dof0) = Eigen::SparseMatrix<double>(2, 2);
    m(dof0, dof1) = Eigen::SparseMatrix<double>(2, 4);
    m(dof1, dof0) = Eigen::SparseMatrix<double>(4, 2);
    m(dof1, dof1) = Eigen::SparseMatrix<double>(4, 4);

    Eigen::SparseMatrix<double> e(6, 6);

    m(dof0, dof0).coeffRef(0, 1) = 1;
    m(dof0, dof0).coeffRef(1, 0) = 10;
    e.coeffRef(0, 1) = 1;
    e.coeffRef(1, 0) = 10;

    m(dof0, dof1).coeffRef(0, 3) = 4;
    e.coeffRef(0, 3 + 2) = 4;

    m(dof1, dof0).coeffRef(2, 1) = 5;
    e.coeffRef(2 + 2, 1) = 5;

    m(dof1, dof1).coeffRef(2, 3) = 6;
    e.coeffRef(2 + 2, 3 + 2) = 6;

    BOOST_CHECK_SMALL((ToEigen(m, {dof0, dof1}) - e).norm(), 1.e-10);
}

BOOST_AUTO_TEST_CASE(DofMatrixZeros)
{
    NuTo::DofType dof0("foo", 1);
    NuTo::DofType dof1("bar", 1);

    NuTo::DofMatrixSparse<double> m;
    m(dof0, dof0) = Eigen::SparseMatrix<double>(2, 8);
    m(dof0, dof1) = Eigen::SparseMatrix<double>(0, 0);
    m(dof1, dof0) = Eigen::SparseMatrix<double>(0, 0);
    m(dof1, dof1) = Eigen::SparseMatrix<double>(0, 10);

    Eigen::SparseMatrix<double> e;
    BOOST_CHECK_NO_THROW(e = ToEigen(m, {dof0, dof1}));
    BOOST_CHECK_EQUAL(e.rows(), 2);
    BOOST_CHECK_EQUAL(e.cols(), 18);
}
