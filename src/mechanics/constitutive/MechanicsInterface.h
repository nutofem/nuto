#pragma once

#include "mechanics/constitutive/EngineeringStrain.h"
#include "mechanics/constitutive/EngineeringStress.h"
#include "mechanics/constitutive/EngineeringTangent.h"
#include "mechanics/cell/CellIds.h"

namespace NuTo
{
namespace Laws
{
template <int TDim>
struct MechanicsInterface
{
    virtual EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double deltaT, CellIds ids) const = 0;
    virtual EngineeringTangent<TDim> Tangent(EngineeringStrain<TDim>, double deltaT, CellIds ids) const = 0;
};
} /* Laws */
} /* NuTo */
