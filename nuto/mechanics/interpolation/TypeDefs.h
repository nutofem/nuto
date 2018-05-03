#pragma once
#include <Eigen/Core>


namespace NuTo
{

constexpr int maxDim = 3;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, NuTo::maxDim, 1> NaturalCoords;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, NuTo::maxDim, 1> GlobalCoords;

} /* NuTo */
