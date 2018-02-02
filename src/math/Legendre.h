#pragma once

#include "base/Exception.h"
#include <vector>

namespace NuTo
{
namespace Math
{

namespace Polynomial
{

//! @brief value of the k-th derivative of legendre polynomial
//! @remark this includes the 0-th derivative, the actual legendre polynomial
//! @param n polynomial order
//! @param x argument of polynomial
//! @param k derivative order
//! @return value of polynomial at x
double Legendre(int n, double x, int k = 0);

//! @brief roots of legendre polynomial
//! @param n polynomial order
//! @return vector of the n-1 roots
std::vector<double> LegendreRoots(int n);

//! @brief roots of first derivative of legendre polynomial
//! @param n polynomial order
//! @return vector of the n-2 roots
std::vector<double> LegendreDerivRoots(int n);

//! @brief computes points and weights for Lobatto quadrature in 1D
//! @param numIPs number of integration points
//! @return pair of quadrature weights and points range [-1,1] including boundary points
std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1DLobatto(int nIps);

//! @brief computes points and weights for Gauss quadrature in 1D
//! @param numIPs number of integration points
//! @return pair of quadrature weights and points range (-1,1)
std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1DGauss(int nIps);

} // namespace Polynomial
} // namespace Math
} // namespace NuTo
