#pragma once

#include <vector>

namespace NuTo
{
namespace Math
{

//! @brief computes points and weights for Lobatto quadrature in 1D
//! @param nIps number of integration points
//! @return pair of quadrature weights and points range [-1,1] including boundary points
std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1DLobatto(int nIps);

//! @brief computes points and weights for Gauss quadrature in 1D
//! @param nIps number of integration points
//! @return pair of quadrature weights and points range (-1,1)
std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1DGauss(int nIps);

} // namespace Math
} // namespace NuTo
