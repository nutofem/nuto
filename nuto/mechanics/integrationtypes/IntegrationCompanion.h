#pragma once

#include <memory>
#include "IntegrationTypeBase.h"
#include "nuto/math/shapes/Shape.h"

namespace NuTo
{

//! Create a new integrationtype (Gauss integration) from a given shape and order.
std::unique_ptr<IntegrationTypeBase> CreateGaussIntegrationType(const Shape& shape, int order);

//! Create a new integrationtype (Lobatto integration) from a given shape and order.
std::unique_ptr<IntegrationTypeBase> CreateLobattoIntegrationType(const Shape& shape, int order);

} // namespace NuTo
