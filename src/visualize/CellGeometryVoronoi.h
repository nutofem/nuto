#pragma once

#include <vector>
#include "mechanics/interpolation/TypeDefs.h"
#include "visualize/VisualizeEnum.h"

namespace NuTo
{
namespace Visualize
{
struct CellGeometryVoronoi
{
    eCellTypes mCellType;
    std::vector<NaturalCoords> mCellCornerCoords;
    //#cells =  mVoronoiCells.size()
    std::vector<std::vector<int>> mVoronoiCells;
};
} /* Visualize */
} /* NuTo */
