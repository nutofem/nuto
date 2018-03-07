#pragma once

#include <memory>
#include "InterpolationSimple.h"
#include "math/shapes/Shape.h"

namespace NuTo
{

//! Create a new serendipity interpolation from a given shape and order.
std::unique_ptr<InterpolationSimple> CreateSerendipityInterpolation(const Shape& shape, int order);

//! Create a new Lobatto interpolation from a given shape and order.
std::unique_ptr<InterpolationSimple> CreateLobattoInterpolation(const Shape& shape, int order);

} // namespace NuTo
