#pragma once

#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIpData.h"

namespace NuTo
{

//! @brief automatically create the lambda
//! [&](cellData, cellIpData) {return integrand.Gradient(cellData, cellIpData, 0); }
template <typename TObject, typename TReturn>
auto Bind(TObject& object, TReturn (TObject::*f)(const NuTo::CellData&, const NuTo::CellIpData&, double))
{
    return [&object, f](const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData) {
        return (object.*f)(cellData, cellIpData, /* deltaT = */ 0.);
    };
}
//! @brief automatically create the lambda
//! [&](cellData, cellIpData) {return integrand.Gradient(cellData, cellIpData); }
template <typename TObject, typename TReturn>
auto Bind(TObject& object, TReturn (TObject::*f)(const NuTo::CellData&, const NuTo::CellIpData&))
{
    return [&object, f](const NuTo::CellData& cellData, const NuTo::CellIpData& cellIpData) {
        return (object.*f)(cellData, cellIpData);
    };
}


} /* NuTo */
