#pragma once

#include <functional>
#include "Eigen/Core"

namespace NuTo
{

class GeometryMeshFem;

namespace UnitMeshFem
{

//! @brief creates a 1 dimensional mesh from (0) -- (1) with numX divisions
//! @param numX number of divisions in x direction
//! @return created mesh
GeometryMeshFem CreateLines(int numX);

//! @brief creates a triangular mesh from (0,0) -- (1,1) with numX and numY divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @return created mesh
GeometryMeshFem CreateTriangles(int numX, int numY);

//! @brief creates a quad mesh from (0,0) -- (1,1) with numX and numY divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @return created mesh
GeometryMeshFem CreateQuads(int numX, int numY);

//! @brief creates a brick mesh from (0,0,0) -- (1,1,1) with numX, numY and numZ divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @param numZ number of divisions in z direction
//! @return created mesh
GeometryMeshFem CreateBricks(int numX, int numY, int numZ);

//! @brief transforms a mesh with a given transformation function f
//! @param oldMesh mesh that is transformed. Call with an xvalue of mesh.
//! @param f transformation function
GeometryMeshFem Transform(GeometryMeshFem&& oldMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f);

} /* UnitGeometryMeshFem */
} /* NuTo */
