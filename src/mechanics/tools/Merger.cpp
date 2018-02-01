#include "mechanics/tools/Merger.h"

using namespace NuTo;

Merger::Merger(MeshFem* rMesh)
    : mMesh(*rMesh)
{
}

Group<NodeSimple>& Merger::Nodes(DofType dof)
{
    if (!mNodes.Has(dof))
        mNodes[dof] = mMesh.NodesTotal(dof);

    return mNodes[dof];
}

void Merger::Merge(GlobalDofVector newValues, std::vector<DofType> dofs)
{
    for (auto dof : dofs)
        for (auto& node : Nodes(dof))
            for (int iDim = 0; iDim < node.GetNumValues(); ++iDim)
            {
                const int dofNumber = node.GetDofNumber(iDim);
                const double dofValue = newValues(dof, dofNumber);
                node.SetValue(iDim, dofValue);
            }
}
