#pragma once

#include "math/LinearInterpolation.h"

namespace NuTo
{
namespace SI
{


//! @brief Density of liquid water at standard atmospheric pressure depending on the temperature.
//!        Values linear interpolated from tables.
//! @param temperature Temperature in Kelvin.
//! @return Density of liquid water.
//! @remark https://en.wikipedia.org/wiki/Density
inline double DensityLiquidWater(double temperature)
{
    const auto densityInterpolation = NuTo::Math::LinearInterpolation({{243.15, 983.854},
                                                                       {253.15, 993.547},
                                                                       {263.15, 998.117},
                                                                       {273.15, 999.8395},
                                                                       {277.15, 999.972},
                                                                       {283.15, 999.7026},
                                                                       {288.15, 999.1026},
                                                                       {293.15, 998.2071},
                                                                       {295.15, 997.7735},
                                                                       {298.15, 997.0479},
                                                                       {303.15, 995.6502},
                                                                       {313.15, 992.2},
                                                                       {333.15, 983.2},
                                                                       {353.15, 971.8},
                                                                       {373.15, 958.4}});
    return densityInterpolation(temperature);
}

} // namespace SI
} // namespace NuTo
