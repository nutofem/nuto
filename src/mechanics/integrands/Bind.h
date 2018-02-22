#pragma once

#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIpData.h"

namespace NuTo
{

//! @brief automatically create the lambda
//! [&](cellIpData) {return integrand.Gradient(cellIpData, 0); }
template <typename TObject, typename TReturn>
auto Bind(TObject& object, TReturn (TObject::*f)(const NuTo::CellIpData&, double))
{
    return [&object, f](const NuTo::CellIpData& cellIpData) { return (object.*f)(cellIpData, /* deltaT = */ 0.); };
}
//! @brief automatically create the lambda
//! [&](cellIpData) {return integrand.Gradient(cellIpData); }
template <typename TObject, typename TReturn>
auto Bind(TObject& object, TReturn (TObject::*f)(const NuTo::CellIpData&))
{
    return [&object, f](const NuTo::CellIpData& cellIpData) { return (object.*f)(cellIpData); };
}


} /* NuTo */
