#pragma once
#include "nuto/mechanics/iga/Nurbs.h"
#include "math/EigenCompanion.h"

namespace NuTo
{
namespace UnitMeshIga
{

//! @brief creates a 1 dimensional mesh from (0) -- (1) with numX divisions (elements)
//! @param numX number of divisions (elements) in x direction
//! @return created mesh = Nurbs geometry
Nurbs<1> CreateLines(int numX)
{
    throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
}

//! @brief creates a quad mesh from (0,0) -- (1,1) with numX and numY divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @return created mesh
Nurbs<2> CreateQuads(int numX, int numY)
{
    throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
}

//! @brief transforms a mesh with a given transformation function f
//! @param oldMesh mesh that is transformed. Call with an xvalue of mesh.
//! @param f transformation function
template <int TDimParameter>
Nurbs<TDimParameter> Transform(Nurbs<TDimParameter>&& oldMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f)
{
    throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
}

} /* UnitMeshFem */
} /* NuTo */
