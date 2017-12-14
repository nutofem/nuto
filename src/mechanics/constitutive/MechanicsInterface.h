#pragma once

#include "mechanics/constitutive/EngineeringStrain.h"
#include "mechanics/constitutive/EngineeringStress.h"

namespace NuTo
{
namespace Laws
{
template <int TDim>
struct MechanicsInterface
{
    using MechanicsTangent = Eigen::Matrix<double, Voigt::Dim(TDim), Voigt::Dim(TDim)>;
    virtual EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double deltaT, int cellId,
                                           int ipId) const = 0;
    virtual MechanicsTangent Tangent(EngineeringStrain<TDim>, double deltaT, int cellId, int ipId) const = 0;
};
} /* Laws */
} /* NuTo */
