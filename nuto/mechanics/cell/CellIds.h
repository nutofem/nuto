#pragma once

namespace NuTo
{

//! named pair for cellId and ipId to make the argument list shorter and avoid accidental mixup of both
struct CellIds
{
    int cellId;
    int ipId;
};
} /* NuTo */
