#include "BoostUnitTest.h"
#include "math/NewtonRaphson.h"
#include <cmath>
#include <iostream>
#include <iomanip>

struct F : NuTo::ResidualDerivative<double, double, double>
{
    double Norm(const double& x) override
    {
        return std::fabs(x);
    }
    
    virtual double R(const double& x) override
    {
        return x * x * x - x + 6; 
    }
    
    virtual double DR(const double& x) override
    {
        return 3. * x * x - 1;
    }

    void Info(unsigned i, const double& x, const double& r) const override
    {
        std::cout << std::setprecision(20);
        std::cout << i << "\t:" << x << "\t" << r << std::endl;
    }
};

struct FInvalid : F 
{
    double R(const double& x) override
    {
        return x * x + 1; 
    }
    
    double DR(const double& x) override
    {
        return 2 * x;
    }
};

struct DoubleSolver
{
    double Solve (double dr, double r) const
    {
        return r / dr;
    }
};

BOOST_AUTO_TEST_CASE(NewtonScalar)
{
    double tolerance = 1.e-10;
    NuTo::NewtonRaphson<F> newton(tolerance, 100);
    F f;
    DoubleSolver solver;
    BOOST_CHECK_CLOSE_FRACTION(newton.Solve(f, 0, solver), -2, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(newton.Solve(f, 0, solver), -2, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalarInvalid)
{
    double tolerance = 1.e-10;
    NuTo::NewtonRaphson<FInvalid> newton(tolerance, 10);
    FInvalid f;
    DoubleSolver solver;
    BOOST_CHECK_THROW(newton.Solve(f, 0, solver), NuTo::NoConvergence);
}

BOOST_AUTO_TEST_CASE(NewtonScalarLineSearch)
{
    double tolerance = 1.e-10;
    NuTo::LineSearchTrue<F> lineSearch(tolerance, 0.01);
    NuTo::NewtonRaphson<F, NuTo::LineSearchTrue<F>> newton(lineSearch, 100);
    F f;
    DoubleSolver solver;
    BOOST_CHECK_CLOSE_FRACTION(newton.Solve(f, 0, solver), -2, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(newton.Solve(f, 0, solver), -2, tolerance);
}
