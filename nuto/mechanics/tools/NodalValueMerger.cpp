#include "nuto/mechanics/tools/NodalValueMerger.h"

using namespace NuTo;

NodalValueMerger::NodalValueMerger(MeshFem* rMesh)
    : mMesh(*rMesh)
{
}

Group<NodeSimple>& NodalValueMerger::Nodes(DofType dof)
{
    if (!mNodes.Has(dof))
        mNodes[dof] = mMesh.NodesTotal(dof);

    return mNodes[dof];
}

void NodalValueMerger::Merge(const GlobalDofVector& newValues, std::vector<DofType> dofs)
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

void NodalValueMerger::Extract(GlobalDofVector* rNewValues, std::vector<DofType> dofs)
{
    for (auto dof : dofs)
        for (auto& node : Nodes(dof))
            for (int iDim = 0; iDim < node.GetNumValues(); ++iDim)
            {
                const int dofNumber = node.GetDofNumber(iDim);
                const double dofValue = node.GetValues()[iDim];
                (*rNewValues)(dof, dofNumber) = dofValue;
            }
}
