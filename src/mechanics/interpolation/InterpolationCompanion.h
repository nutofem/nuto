#pragma once

#include <memory>
#include "InterpolationSimple.h"
#include "math/shapes/Shape.h"

namespace NuTo
{

//! Create a new interpolation from a given shape and order.
std::unique_ptr<InterpolationSimple> CreateInterpolation(const Shape& shape, int order);

} // namespace NuTo
