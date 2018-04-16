#pragma once
#include "nuto/mechanics/dofs/GlobalDofVector.h"
#include "nuto/mechanics/mesh/MeshFem.h"

namespace NuTo
{
//! Performs our good old "NodeMerge" and should be replaced by any solution from issue #141 PDE nodal values
class NodalValueMerger
{
public:
    NodalValueMerger(MeshFem* rMesh);

    //! Performs "NodeMerge", writes values from the global solution vector to the nodes
    //! @param newValues complete vector (J/K) of new dof values
    //! @param dofs dof types to merge
    //! @param instance id of the dof instance
    void Merge(const GlobalDofVector& newValues, std::vector<DofType> dofs, int instance = 0);

    //! Performs "NodeExtract", writes values from the global solution vector to the nodes
    //! @param rNewValues complete vector (J/K) of new dof values
    //! @param dofs dof types to extract
    //! @param instance id of the dof instance
    void Extract(GlobalDofVector* rNewValues, std::vector<DofType> dofs, int instance = 0);

    //! node group memoizer
    //! @param dof dof type
    //! @return a memoized node group
    Group<DofNode>& Nodes(DofType dof);

private:
    MeshFem& mMesh;
    DofContainer<Group<DofNode>> mNodes;
};
} /* NuTo */
