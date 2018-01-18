#pragma once

#include "mechanics/mesh/MeshFem.h"

namespace NuTo
{
//! @brief add a new _layer_ of nodes for `dofType`
//! @param rMesh fem mesh, return argument with r and wierd pointer syntax to make it clear
//! @param dofType dof type
//! @param interpolation interpolation type for the dof type
void AddDofInterpolation(NuTo::MeshFem* rMesh, DofType dofType, const InterpolationSimple& interpolation);

//! @brief change coordinate interpolation
//! @param rMesh fem mesh, return argument with r and wierd pointer syntax to make it clear
//! @param interpolation interpolation type for the dof type
//! @remark New nodes are created, the old ones are still there, the same holds for the interpolation
void ChangeCoordinateInterpolation(NuTo::MeshFem* rMesh, const InterpolationSimple& interpolation);

} /* NuTo */
