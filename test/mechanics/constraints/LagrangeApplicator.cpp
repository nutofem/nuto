#include "BoostUnitTest.h"
#include "Eigen/Sparse"
#include "Eigen/SVD"
#include <iostream>
#include <utility>

#include "nuto/mechanics/constraints/LagrangeApplicator.h"
#include "nuto/math/EigenSparseSolve.h"

using namespace NuTo;

Eigen::SparseMatrix<double> TrussSystem(int n)
{
    Eigen::SparseMatrix<double> K;
    K.resize(n+1, n+1);
    K.insert(0, 0) = 1.0;
    for (int i = 1; i < n+1; i++)
    {
        K.insert(i, i) = 2.0;
        K.insert(i, i-1) = -1.0;
        K.insert(i-1, i) = -1.0;
    }
    K.insert(n, n) = 1.0;
    K.makeCompressed();
    return K;
}

struct ConstrainedSystem
{
    Math::Equation equation;
    Eigen::VectorXd f;
    Eigen::VectorXd expected_result;
    Eigen::SparseMatrix<double> K = TrussSystem(2);
    
    ConstrainedSystem()
    {
        Math::Term u3(2, 1.0);
        equation.Terms().push_back(u3);
        equation.Value() = 1.0;

        f.resize(3);
        f << 1.0, 2.0, 0.0;

        expected_result.resize(3);
        expected_result << 5.0, 4.0, 1.0;
    }
};

struct PeriodicSystem
{
    Math::Equation periodicEquation;
    Math::Equation centerEquation;
    Eigen::VectorXd f;
    Eigen::VectorXd expected_result;
    Eigen::SparseMatrix<double> K = TrussSystem(4);
    PeriodicSystem()
    {

        f.resize(5);
        f.setZero();

        Math::Term u1(0, -1.0);
        Math::Term u5(4, 1.0);
        periodicEquation.Terms().push_back(u1);
        periodicEquation.Terms().push_back(u5);
        periodicEquation.Value() = 1.0;

        Math::Term u3(2, 1.0);
        centerEquation.Terms().push_back(u3);
        centerEquation.Value() = 0.0;

        expected_result.resize(5);
        expected_result << -0.5, -0.25, 0.0, 0.25, 0.5;
    }
};


BOOST_FIXTURE_TEST_CASE(ApplyUsingLagrange, ConstrainedSystem)
{
    LagrangeApplicator applicator({equation});

    auto K_mod = applicator.Apply(K);
    auto f_mod = applicator.Apply(f);
    auto u_mod = EigenSparseSolve(K_mod, f_mod, "EigenSparseLU");

    auto u = applicator.ConvertBack(u_mod);

    BoostUnitTest::CheckEigenMatrix(u, expected_result);
}

BOOST_FIXTURE_TEST_CASE(PeriodicUsingLagrange, PeriodicSystem)
{
    LagrangeApplicator applicator({periodicEquation, centerEquation});

    auto K_mod = applicator.Apply(K);
    auto f_mod = applicator.Apply(f);
    auto u_mod = EigenSparseSolve(K_mod, f_mod, "EigenSparseLU");

    auto u = applicator.ConvertBack(u_mod);

    BoostUnitTest::CheckEigenMatrix(u, expected_result);
}
