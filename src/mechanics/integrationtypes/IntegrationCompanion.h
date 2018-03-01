#pragma once

#include <memory>
#include "IntegrationTypeBase.h"
#include "math/shapes/Shape.h"

namespace NuTo
{

//! Create a new integrationtype from a given shape and order.
std::unique_ptr<IntegrationTypeBase> CreateIntegrationType(const Shape& shape, int order);

} // namespace NuTo
