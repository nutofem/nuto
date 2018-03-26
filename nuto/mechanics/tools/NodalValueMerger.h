#pragma once
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/mesh/MeshFem.h"

namespace NuTo
{
//! Performs our good old "NodeMerge" and should be replaced by any solution from issue #141 PDE nodal values
class NodalValueMerger
{
public:
    NodalValueMerger(MeshFem* rMesh);

    //! Performs "NodeMerge", writes values from the solution vector to the nodes
    //! @param new dof values
    //! @param dofs dof types to merge
    void Merge(const DofVector<double>& newValues, std::vector<DofType> dofs);

    //! Performs "NodeExtract", writes values from the solution vector to the nodes
    //! @param rNewValues
    //! @param dofs dof types to extract
    void Extract(DofVector<double>* rNewValues, std::vector<DofType> dofs);

    //! node group memoizer
    //! @param dof dof type
    //! @return a memoized node group
    Group<NodeSimple>& Nodes(DofType dof);

private:
    MeshFem& mMesh;
    DofContainer<Group<NodeSimple>> mNodes;
};
} /* NuTo */
