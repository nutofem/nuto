#pragma once

#include <vector>
#include "mechanics/interpolation/TypeDefs.h"
#include "visualize/VisualizeEnum.h"

namespace NuTo
{
namespace Visualize
{
struct CellGeometry
{
    eCellTypes mCellType;
    std::vector<NaturalCoords> mCornerCoords;
    std::vector<NaturalCoords> mIntegrationPointCoords;
};
} /* Visualize */
} /* NuTo */
