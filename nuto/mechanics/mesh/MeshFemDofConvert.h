#pragma once

#include "nuto/mechanics/mesh/MeshFem.h"
#include <boost/optional.hpp>

namespace NuTo
{
//! @brief add a new _layer_ of nodes for `dofType`
//! @param rMesh fem mesh, return argument with r and weird pointer syntax to make it clear
//! @param dofType dof type
//! @param optionalInterpolation interpolation type for the dof type. If you do _not_ provide a argument here (or
//! boost::none) each new dof element will use the underlying coordinate interpolation resulting in an isoparametric
//! element. This will also work for mixed meshes consisting of multiple interpolation types
void AddDofInterpolation(NuTo::MeshFem* rMesh, DofType dofType,
                         boost::optional<const InterpolationSimple&> optionalInterpolation = boost::none);

//! @brief add a new _layer_ of nodes for `dofType` but only for elements in given group
void AddDofInterpolation(NuTo::MeshFem* rMesh, DofType dofType, Group<ElementCollectionFem> elements,
                         boost::optional<const InterpolationSimple&> optionalInterpolation = boost::none);

} /* NuTo */
