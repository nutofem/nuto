#pragma once

#include "visualize/VoronoiHandler.h"

namespace NuTo
{
namespace Visualize
{

//! Creates a voronoi geometry for equidistant lines in 1D
//! @param numCellsPerDirection Number of cells "per Direction" i.e. 3 --> 3 voronoi cells
VoronoiGeometry VoronoiGeometryLine(int numCellsPerDirection);

//! Creates a voronoi geometry for equidistant quads in 2D
//! @param numCellsPerDirection Number of cells "per Direction" i.e. 3 --> 9 num voronoi cells
VoronoiGeometry VoronoiGeometryQuad(int numCellsPerDirection);

//! Creates a voronoi geometry for equidistant bricks in 3D
//! @param numCellsPerDirection Number of cells "per Direction" i.e. 3 --> 27 num voronoi cells
VoronoiGeometry VoronoiGeometryBrick(int numCellsPerDirection);

//! Creates a voronoi geometry for a triangle with 3 integration points at 1/6* ( 1,1; 4,1; 1,4 )
VoronoiGeometry VoronoiGeometryTriangle3Ip();

} // namespace Visualize
} // namespace NuTo
