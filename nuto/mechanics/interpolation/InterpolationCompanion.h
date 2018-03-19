#pragma once

#include <memory>
#include "InterpolationSimple.h"
#include "nuto/math/shapes/Shape.h"

namespace NuTo
{

//! Create a new Lagrange interpolation from a given shape and order.
std::unique_ptr<InterpolationSimple> CreateLagrangeInterpolation(const Shape& shape, int order);

//! Create a new Lobatto interpolation from a given shape and order.
std::unique_ptr<InterpolationSimple> CreateLobattoInterpolation(const Shape& shape, int order);

} // namespace NuTo
