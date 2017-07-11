#pragma once
#include <string>
#include <map>
#include "visualize/VisualizeEnum.h"

namespace NuTo
{

std::string GetComponentName(eVisualizeWhat component);

eVisualizeWhat GetComponentEnum(std::string component);

} /* NuTo */
