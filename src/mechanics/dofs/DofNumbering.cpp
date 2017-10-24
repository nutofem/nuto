#include "mechanics/dofs/DofNumbering.h"

using namespace NuTo;

//! @brief build dof numbering, starting at 0, for all `nodes` regardless of constraints
//! @return total number of dofs in `nodes`
int InitialUnconstrainedNumbering(const Groups::Group<NodeSimple>& nodes)
{
    int dofNumber = 0;
    for (auto& node : nodes)
        for (int iComponent = 0; iComponent < node.GetNumValues(); ++iComponent)
            node.SetDofNumber(iComponent, dofNumber++);
    return dofNumber;
}

std::vector<bool> FindConstrainedDofs(const ConstraintPde::Constraints& constraints, DofType dof, int numDofs)
{
    std::vector<bool> isConstrained(numDofs, false);
    for (int iEquation = 0; iEquation < constraints.GetNumEquations(dof); ++iEquation)
    {
        int dependentDofNumber = constraints.GetEquation(dof, iEquation).GetDependentDofNumber();
        isConstrained[dependentDofNumber] = true;
    }
    return isConstrained;
}

DofNumbering::DofInfo DofNumbering::Build(const Groups::Group<NodeSimple>& dofNodes, DofType dof,
                                          const ConstraintPde::Constraints& constraints)
{
    int numDependentDofs = constraints.GetNumEquations(dof);

    const int numDofs = InitialUnconstrainedNumbering(dofNodes);

    std::vector<bool> isConstrained = FindConstrainedDofs(constraints, dof, numDofs);

    int countIndependentDofs = 0;
    int countDependentDofs = numDofs - numDependentDofs;
    for (auto& node : dofNodes)
        for (int iComponent = 0; iComponent < node.GetNumValues(); ++iComponent)
        {
            int dofNumber = node.GetDofNumber(iComponent);
            if (isConstrained[dofNumber])
                node.SetDofNumber(iComponent, countDependentDofs++);
            else
                node.SetDofNumber(iComponent, countIndependentDofs++);
        }

    DofInfo dofInfo;
    dofInfo.numDependentDofs[dof] = numDependentDofs;
    dofInfo.numIndependentDofs[dof] = numDofs - numDependentDofs;
    return dofInfo;
}
