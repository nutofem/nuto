#include "BoostUnitTest.h"
#include <Eigen/Eigenvalues>
#include <iostream> // ArpackSupport writes to std::cout without including <iostream> (bug?). So we have to include it.
#include <Eigen/Sparse>

#include <unsupported/Eigen/ArpackSupport>

Eigen::SparseMatrix<double> GetA()
{
    Eigen::SparseMatrix<double> a(8, 8);
    a.insert(0, 0) = 8.;
    a.insert(1, 1) = 2.;
    a.insert(2, 2) = 4.;
    a.insert(3, 3) = 7.;
    a.insert(4, 4) = 10.;
    a.insert(5, 5) = 6.;
    a.insert(6, 6) = 5.;
    a.insert(7, 7) = 1.;
    a.insert(0, 1) = 1.;
    a.insert(1, 0) = 1.;
    a.insert(0, 3) = -1.;
    a.insert(3, 0) = -1.;
    a.insert(2, 4) = 2.;
    a.insert(4, 2) = 2.;
    a.insert(6, 1) = -3.;
    a.insert(1, 6) = -3.;
    return a;
}

Eigen::SparseMatrix<double> GetM()
{
    Eigen::SparseMatrix<double> m(8, 8);
    m.insert(0, 0) = 1.;
    m.insert(1, 1) = 1.;
    m.insert(2, 2) = 1.;
    m.insert(3, 3) = 1.;
    m.insert(4, 4) = 1.;
    m.insert(5, 5) = 1.;
    m.insert(6, 6) = 1.;
    m.insert(7, 7) = 1.;

    m.insert(0, 1) = 0.5;
    m.insert(1, 0) = 0.5;
    return m;
}

Eigen::VectorXd EVs(Eigen::MatrixXd A, Eigen::MatrixXd M)
{
    return Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>(A, M).eigenvalues();
}

Eigen::VectorXd EVs(Eigen::MatrixXd A)
{
    return Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(A).eigenvalues();
}

BOOST_AUTO_TEST_CASE(ArpackUsage)
{
    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> arpack;
    Eigen::SparseMatrix<double> A = GetA();

    // find smallest EV: - smallest magnitude ("SM") or smallest amplitude ("SA")
    BOOST_CHECK_CLOSE(arpack.compute(A, 1, "SM", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A).minCoeff(), 1.e-10);
    BOOST_CHECK_CLOSE(arpack.compute(A, 1, "SA", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A).minCoeff(), 1.e-10);

    // find largest EV: - largest magnitude ("LM") or largest amplitude ("LA")
    BOOST_CHECK_CLOSE(arpack.compute(A, 1, "LM", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A).maxCoeff(), 1.e-10);
    BOOST_CHECK_CLOSE(arpack.compute(A, 1, "LA", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A).maxCoeff(), 1.e-10);
}

BOOST_AUTO_TEST_CASE(ArpackUsageGeneral)
{
    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> arpack;
    Eigen::SparseMatrix<double> A = GetA();
    Eigen::SparseMatrix<double> M = GetM();

    // find smallest EV: - smallest magnitude ("SM") or smallest amplitude ("SA")
    BOOST_CHECK_CLOSE(arpack.compute(A, M, 1, "SM", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A, M).minCoeff(),
                      1.e-10);
    BOOST_CHECK_CLOSE(arpack.compute(A, M, 1, "SA", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A, M).minCoeff(),
                      1.e-10);

    // find largest EV: - largest magnitude ("LM") or largest amplitude ("LA")
    BOOST_CHECK_CLOSE(arpack.compute(A, M, 1, "LM", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A, M).maxCoeff(),
                      1.e-10);
    BOOST_CHECK_CLOSE(arpack.compute(A, M, 1, "LA", Eigen::EigenvaluesOnly).eigenvalues()[0], EVs(A, M).maxCoeff(),
                      1.e-10);
}
