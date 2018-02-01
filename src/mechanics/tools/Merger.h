#pragma once
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/mesh/MeshFem.h"

namespace NuTo
{
//! Performs our good old "NodeMerge" and should be replaced by any solution from issue #141 PDE nodal values
class Merger
{
public:
    Merger(MeshFem* rMesh);

    //! Performs "NodeMerge", writes values from the global solution vector to the nodes
    //! @param newValues complete vector (J/K) of new dof values
    //! @param dofs dof types to merge
    void Merge(GlobalDofVector newValues, std::vector<DofType> dofs);

private:
    MeshFem& mMesh;
    DofContainer<Group<NodeSimple>> mNodes;

    //! node group memoizer
    //! @param dof dof type
    //! @return a memoized node group
    Group<NodeSimple>& Nodes(DofType dof);
};
} /* NuTo */
