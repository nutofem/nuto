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
//! @param rMesh mesh that is modified (return argument! pointer syntax to make it clear)
//! @param f transformation function
void Transform(MeshFem* rMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f);

} /* UnitMeshFem */
} /* NuTo */
