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
    virtual EngineeringStressPDE<TDim> Stress(EngineeringStrainPDE<TDim> strain, double deltaT) const = 0;
    virtual MechanicsTangent Tangent(EngineeringStrainPDE<TDim>, double deltaT) const = 0;
};
} /* Laws */
} /* NuTo */
