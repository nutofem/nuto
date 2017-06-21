#include "BoostUnitTest.h"
#include "math/NewtonRaphson.h"
#include <cmath>
//#include <iostream>
//#include <iomanip>

constexpr double tolerance = 1.e-10;

auto Problem()
{
    auto R = [](double x) { return x * x * x - x + 6; };
    auto DR = [](double x) { return 3. * x * x - 1; };
    auto Norm = [](double x) { return std::abs(x); };
    return NuTo::DefineNonlinearProblem(R, DR, Norm, tolerance);
}

auto InvalidProblem()
{
    auto R = [](double x) { return x * x + 1; };
    auto DR = [](double x) { return 2. * x; };
    auto Norm = [](double x) { return std::abs(x); };
    return NuTo::DefineNonlinearProblem(R, DR, Norm, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalar)
{
    auto result = NuTo::Newton(Problem(), 0., NuTo::DoubleSolver(), 100);
    BOOST_CHECK_CLOSE_FRACTION(result, -2, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalarLineSearch)
{
    auto result = NuTo::Newton(Problem(), 0., NuTo::DoubleSolver(), 20, NuTo::LineSearch());
    BOOST_CHECK_CLOSE_FRACTION(result, -2, tolerance);
}

BOOST_AUTO_TEST_CASE(NewtonScalarInvalid)
{
    BOOST_CHECK_THROW(NuTo::Newton(InvalidProblem(), 0., NuTo::DoubleSolver(), 100), NuTo::NoConvergence);
}

BOOST_AUTO_TEST_CASE(NewtonScalarLineSearchInvalid)
{
    BOOST_CHECK_THROW(NuTo::Newton(InvalidProblem(), 0., NuTo::DoubleSolver(), 20, NuTo::LineSearch()),
                      NuTo::NoConvergence);
}
