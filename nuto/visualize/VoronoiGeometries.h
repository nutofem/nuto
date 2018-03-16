#pragma once

#include "nuto/visualize/VoronoiHandler.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
namespace Visualize
{
enum Spacing
{
    EQUIDISTANT,
    LOBATTO,
    GAUSS
};

//! Creates a voronoi geometry for equidistant lines in 1D
//! @param numCellsPerDirection Number of cells "per Direction" i.e. 3 --> 3 voronoi cells
VoronoiGeometry VoronoiGeometryLine(int numCellsPerDirection, Spacing s = EQUIDISTANT);

//! Creates a voronoi geometry for equidistant quads in 2D
//! @param numCellsPerDirection Number of cells "per Direction" i.e. 3 --> 9 num voronoi cells
VoronoiGeometry VoronoiGeometryQuad(int numCellsPerDirection, Spacing s = EQUIDISTANT);

//! Creates a voronoi geometry for equidistant bricks in 3D
//! @param numCellsPerDirection Number of cells "per Direction" i.e. 3 --> 27 num voronoi cells
VoronoiGeometry VoronoiGeometryBrick(int numCellsPerDirection, Spacing s = EQUIDISTANT);

//! Creates a voronoi geometry for an arbitrary triangle integration type
VoronoiGeometry VoronoiGeometryTriangle(const IntegrationTypeBase& integrationType);

} // namespace Visualize
} // namespace NuTo
