#pragma once

#include "mechanics/constitutive/EngineeringStrainPDE.h"
#include "mechanics/constitutive/EngineeringStressPDE.h"

namespace NuTo
{
namespace Laws
{
template <int TDim>
struct MechanicsInterface
{
    using MechanicsTangent = Eigen::Matrix<double, Voigt::Dim(TDim), Voigt::Dim(TDim)>;
    virtual EngineeringStressPDE<TDim> Stress(EngineeringStrainPDE<TDim> strain, double deltaT, int cellId,
                                              int ipNum) const = 0;
    virtual MechanicsTangent Tangent(EngineeringStrainPDE<TDim>, double deltaT, int cellId, int ipId) const = 0;
};
} /* Laws */
} /* NuTo */
