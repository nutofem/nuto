#include "BoostUnitTest.h"
#include <Eigen/Eigenvalues>
#include <iostream>
#include <Eigen/Sparse>

#include <unsupported/Eigen/ArpackSupport>

struct ArpackTestFixture
{
    ArpackTestFixture()
    {
        BOOST_TEST_MESSAGE("setup fixture");
        auto A_Full = GetA();
        auto M_Full = GetM();

        mEvs = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(A_Full).eigenvalues();
        mEvsGeneral = Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>(A_Full, M_Full).eigenvalues();
    }

    static Eigen::SparseMatrix<double> GetA()
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

    static Eigen::SparseMatrix<double> GetM()
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

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> mArpack;

    Eigen::VectorXd mEvs;
    Eigen::VectorXd mEvsGeneral;
};

BOOST_FIXTURE_TEST_CASE(Info, ArpackTestFixture)
{
    BOOST_TEST_MESSAGE("FullMatrix\n" << GetA());
    BOOST_TEST_MESSAGE("EVs of standard eigenvalue problem.\n" << mEvs.transpose());
    BOOST_TEST_MESSAGE("EVs of generalized eigenvalue problem.\n" << mEvsGeneral.transpose());
}

BOOST_FIXTURE_TEST_CASE(LargestAmplitude, ArpackTestFixture)
{
    auto solution = mArpack.compute(GetA(), 3, "LA", Eigen::EigenvaluesOnly);
    BoostUnitTest::CheckEigenMatrix(solution.eigenvalues(), mEvs.segment(5, 3));
}

BOOST_FIXTURE_TEST_CASE(SmallestMagnitude, ArpackTestFixture)
{
    auto solution = mArpack.compute(GetA(), 3, "SM", Eigen::EigenvaluesOnly);
    BoostUnitTest::CheckEigenMatrix(solution.eigenvalues(), mEvs.segment(0, 3));
}

BOOST_FIXTURE_TEST_CASE(LargestAmplitudeGeneral, ArpackTestFixture)
{
    auto solution = mArpack.compute(GetA(), GetM(), 3, "LA", Eigen::EigenvaluesOnly);
    BoostUnitTest::CheckEigenMatrix(solution.eigenvalues(), mEvsGeneral.segment(5, 3));
}

BOOST_FIXTURE_TEST_CASE(SmallestMagnitudeGeneral, ArpackTestFixture)
{
    auto solution = mArpack.compute(GetA(), GetM(), 3, "SM", Eigen::EigenvaluesOnly);
    BoostUnitTest::CheckEigenMatrix(solution.eigenvalues(), mEvsGeneral.segment(0, 3));
}
