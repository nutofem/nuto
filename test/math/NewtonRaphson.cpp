#include "BoostUnitTest.h"
#include "math/NewtonRaphson.h"
#include <cmath>
#include <iostream>
#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseMatrixCSRGeneral.h"
#include <complex>

constexpr double tolerance = 1.e-10;
using namespace NuTo::NewtonRaphson;

/* ##################################################
 * ##              SCALAR REAL TESTS               ##
 * ################################################## */

auto ValidProblem()
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };
    auto Info = [](int i, double x, double r) { std::cout << i << '\t' << x << '\t' << r << '\n'; };
    return DefineProblem(R, DR, Norm, tolerance, Info);
}

template <typename T>
auto InvalidProblem()
{
    auto R = [](T x) { return x * x + 1.; };
    auto DR = [](T x) { return 2. * x; };
    auto Norm = [](T x) { return std::abs(x); };
    return DefineProblem(R, DR, Norm, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalar)
{
    auto result = Solve(ValidProblem(), 0., DoubleSolver(), 100);
    BOOST_CHECK_CLOSE_FRACTION(result, -2, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalarLineSearch)
{
    int numIterations = 0;
    auto result = Solve(ValidProblem(), 0., DoubleSolver(), 20, LineSearch(), &numIterations);
    BOOST_CHECK_CLOSE_FRACTION(result, -2, tolerance);
    BOOST_CHECK(numIterations > 0);
}

BOOST_AUTO_TEST_CASE(NewtonScalarInvalid)
{
    int numIterations = 0;
    BOOST_CHECK_THROW(Solve(InvalidProblem<double>(), 0., DoubleSolver(), 100, NoLineSearch(), &numIterations), NoConvergence);
    BOOST_CHECK_EQUAL(numIterations, 100);
}

BOOST_AUTO_TEST_CASE(NewtonScalarLineSearchInvalid)
{
    BOOST_CHECK_THROW(Solve(InvalidProblem<double>(), 0., DoubleSolver(), 20, LineSearch()), NoConvergence);
}

/* ##################################################
 * ##            SCALAR COMPLEX TESTS              ##
 * ################################################## */

struct ComplexSolver
{
    std::complex<double> Solve(const std::complex<double>& DR, const std::complex<double>& R) const
    {
        return R / DR;
    }
};

BOOST_AUTO_TEST_CASE(NewtonScalarComplex)
{
    auto root1 = Solve(InvalidProblem<std::complex<double>>(), std::complex<double>(0, .1), ComplexSolver(), 100);
    BOOST_CHECK_SMALL(std::abs(root1 - std::complex<double>(0.,1.)), tolerance);
    
    auto root2 = Solve(InvalidProblem<std::complex<double>>(), std::complex<double>(0, -.1), ComplexSolver(), 100);
    BOOST_CHECK_SMALL(std::abs(root2 - std::complex<double>(0.,-1.)), tolerance);
}

/* ##################################################
 * ##              VECTOR REAL TESTS               ##
 * ################################################## */

auto ValidMatrixProblem()
{
    auto R = [](Eigen::VectorXd x) {
        Eigen::VectorXd r(2);
        r[0] = x[0] * x[0] * x[0] - x[0] + 6;
        r[1] = x[1] - 1;
        return r;
    };
    auto DR = [](Eigen::VectorXd x) {
        NuTo::SparseMatrixCSRGeneral<double> m(2, 2);
        m.AddValue(0, 0, 3 * x[0] * x[0] - 1);
        m.AddValue(1, 1, 1);
        return m;
    };
    auto Norm = [](Eigen::VectorXd x) { return x.norm(); };
    return DefineProblem(R, DR, Norm, tolerance);
}

struct MumpsWrapper
{
    Eigen::VectorXd Solve(NuTo::SparseMatrixCSRGeneral<double>& DR, const Eigen::VectorXd& R)
    {
        Eigen::VectorXd v;
        DR.SetOneBasedIndexing();
        mSolver.Solve(DR, R, v);
        return v;
    }
    NuTo::SparseDirectSolverMUMPS mSolver;
};

BOOST_AUTO_TEST_CASE(NewtonSparse)
{
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(2);
    auto result = Solve(ValidMatrixProblem(), x0, MumpsWrapper(), 100);
    BoostUnitTest::CheckVector(result, std::vector<double>{-2., 1.}, 2);
}

BOOST_AUTO_TEST_CASE(NewtonSparseLineSearch)
{
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(2);
    auto result = Solve(ValidMatrixProblem(), x0, MumpsWrapper(), 20, LineSearch());
    BoostUnitTest::CheckVector(result, std::vector<double>{-2., 1.}, 2);
}
