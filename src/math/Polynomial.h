#pragma once

#include <vector>

namespace NuTo {
namespace Math {

namespace Polynomial {

//! @brief value of legendre polynomial
//! @param n polynomial order
//! @param x argument of polynomial
//! @return value of polynomial at x
double Legendre(int n, double x);

//! @brief value of the k-th derivative of legendre polynomial
//! @param n polynomial order
//! @param x argument of polynomial
//! @param k derivative order
//! @return value of polynomial at x
double LegendreDeriv(int n, double x, int k);

//! @brief roots of legendre polynomial
//! @param n polynomial order
//! @return vector of the n-1 roots
std::vector<double> LegendreRoots(int n);

//! @brief roots of first derivative of legendre polynomial
//! @param n polynomial order
//! @return vector of the n-2 roots
std::vector<double> LegendreDerivRoots(int n);

} // namespace Polynomial
} // namespace Math
} // namespace NuTo
