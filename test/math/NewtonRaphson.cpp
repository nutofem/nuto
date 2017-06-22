#include "BoostUnitTest.h"
#include "math/NewtonRaphson.h"
#include <cmath>
#include <iostream>

constexpr double tolerance = 1.e-10;
using namespace NuTo::NewtonRaphson;
auto ValidProblem()
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };
    auto Info = [](int i, double x, double r) { std::cout << i << '\t' << x << '\t' << r << '\n'; };
    return DefineProblem(R, DR, Norm, tolerance, Info);
}

auto InvalidProblem()
{
    auto R = [](double x) { return x * x + 1; };
    auto DR = [](double x) { return 2. * x; };
    auto Norm = [](double x) { return std::abs(x); };
    return DefineProblem(R, DR, Norm, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalar)
{
    auto result = Solve(ValidProblem(), 0., DoubleSolver(), 100);
    BOOST_CHECK_CLOSE_FRACTION(result, -2, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalarLineSearch)
{
    auto result = Solve(ValidProblem(), 0., DoubleSolver(), 20, LineSearch());
    BOOST_CHECK_CLOSE_FRACTION(result, -2, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalarInvalid)
{
    BOOST_CHECK_THROW(Solve(InvalidProblem(), 0., DoubleSolver(), 100), NoConvergence);
}

BOOST_AUTO_TEST_CASE(NewtonScalarLineSearchInvalid)
{
    BOOST_CHECK_THROW(Solve(InvalidProblem(), 0., DoubleSolver(), 20, LineSearch()), NoConvergence);
}
