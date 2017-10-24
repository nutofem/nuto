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
    std::vector<bool> isDependent(numDofs, false);
    for (int iEquation = 0; iEquation < constraints.GetNumEquations(dof); ++iEquation)
    {
        int dependentDofNumber = constraints.GetEquation(dof, iEquation).GetDependentDofNumber();
        if (isDependent[dependentDofNumber])
            throw NuTo::Exception(
                    __PRETTY_FUNCTION__,
                    "Dof number " + std::to_string(dependentDofNumber) +
                            " already constrained. You are only allowed to have a dependent dof in one equation.");
        isDependent[dependentDofNumber] = true;
    }

    // go through all terms that are not dependent and check if they contain dependent dofs
    for (int iEquation = 0; iEquation < constraints.GetNumEquations(dof); ++iEquation)
    {
        for (int iTerm = 1; iTerm < constraints.GetEquation(dof, iEquation).GetTerms().size(); ++iTerm)
        {
            const ConstraintPde::Term& term = constraints.GetEquation(dof, iEquation).GetTerms()[iTerm];
            int termDofNumber = term.GetConstrainedDofNumber();

            if (isDependent[termDofNumber])
                throw NuTo::Exception(__PRETTY_FUNCTION__,
                                      "Term " + std::to_string(iTerm) + " of Equation " + std::to_string(iEquation) +
                                              " contains dof number " + std::to_string(termDofNumber) +
                                              " that is a dependent dof number. This is not allowed!");
        }
    }

    return isDependent;
}

DofNumbering::DofInfo DofNumbering::Build(const Groups::Group<NodeSimple>& dofNodes, DofType dof,
                                          const ConstraintPde::Constraints& constraints)
{
    int numDependentDofs = constraints.GetNumEquations(dof);

    const int numDofs = InitialUnconstrainedNumbering(dofNodes);

    std::vector<bool> isDependent = FindConstrainedDofs(constraints, dof, numDofs);

    int countIndependentDofs = 0;
    int countDependentDofs = numDofs - numDependentDofs;
    for (auto& node : dofNodes)
        for (int iComponent = 0; iComponent < node.GetNumValues(); ++iComponent)
        {
            int dofNumber = node.GetDofNumber(iComponent);
            if (isDependent[dofNumber])
                node.SetDofNumber(iComponent, countDependentDofs++);
            else
                node.SetDofNumber(iComponent, countIndependentDofs++);
        }

    DofInfo dofInfo;
    dofInfo.numDependentDofs[dof] = numDependentDofs;
    dofInfo.numIndependentDofs[dof] = numDofs - numDependentDofs;
    return dofInfo;
}
