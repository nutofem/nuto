#pragma once
#include "mechanics/mesh/MeshFem.h"

namespace NuTo
{
namespace UnitMeshFem
{

//! @brief creates a triangular mesh from (0,0) -- (1,1) with numX and numY divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @return created mesh
MeshFem CreateTriangles(int numX, int numY);

//! @brief creates a quad mesh from (0,0) -- (1,1) with numX and numY divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @return created mesh
MeshFem CreateQuads(int numX, int numY);

//! @brief transforms a mesh with a given transformation function f
//! @param oldMesh mesh that is transformed. Call with an xvalue of mesh.
//! @param f transformation function
MeshFem Transform(MeshFem&& oldMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f);

} /* UnitMeshFem */
} /* NuTo */
