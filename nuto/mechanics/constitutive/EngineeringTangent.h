#pragma once

#include <eigen3/Eigen/Core>
#include "nuto/mechanics/constitutive/Voigt.h"

namespace NuTo
{
template <int TDim>
using EngineeringTangent = Eigen::Matrix<double, Voigt::Dim(TDim), Voigt::Dim(TDim)>;
} /* NuTo */
