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
    //! @param newValues new dof values
    //! @param dofs dof types to merge
    //! @param instance id of the dof instance
    void Merge(const DofVector<double>& newValues, std::vector<DofType> dofs, int instance = 0);

    //! Performs "NodeExtract", writes values from the solution vector to the nodes
    //! @param rNewValues vector to fill the new values with
    //! @param dofs dof types to extract
    //! @param instance id of the dof instance
    void Extract(DofVector<double>* rNewValues, std::vector<DofType> dofs, int instance = 0);

    //! node group memoizer
    //! @param dof dof type
    //! @return a memoized node group
    Group<NodeSimple>& Nodes(DofType dof);

private:
    MeshFem& mMesh;
    DofContainer<Group<NodeSimple>> mNodes;
};
} /* NuTo */
