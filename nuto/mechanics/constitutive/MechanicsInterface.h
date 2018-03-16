#pragma once

#include "nuto/mechanics/constitutive/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/EngineeringStress.h"
#include "nuto/mechanics/cell/CellIds.h"

namespace NuTo
{
namespace Laws
{
template <int TDim>
struct MechanicsInterface
{
    using MechanicsTangent = Eigen::Matrix<double, Voigt::Dim(TDim), Voigt::Dim(TDim)>;
    virtual EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double deltaT, CellIds ids) const = 0;
    virtual MechanicsTangent Tangent(EngineeringStrain<TDim>, double deltaT, CellIds ids) const = 0;
};
} /* Laws */
} /* NuTo */
