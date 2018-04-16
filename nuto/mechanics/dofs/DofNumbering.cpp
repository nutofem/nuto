#include "nuto/mechanics/dofs/DofNumbering.h"

using namespace NuTo;

//! @brief build dof numbering, starting at 0, for all `nodes` regardless of constraints
//! @return total number of dofs in `nodes`
int InitialUnconstrainedNumbering(const Group<DofNode>& nodes)
{
    int dofNumber = 0;
    for (auto& node : nodes)
        for (int iComponent = 0; iComponent < node.GetNumValues(); ++iComponent)
            node.SetDofNumber(iComponent, dofNumber++);
    return dofNumber;
}

constexpr int INDEPENDENT = -1;

//! @brief Find the status (INDEPENDENT or dependent) for all dof numbers.
//! @return The status vector contains two information. Value INDEPENDENT means, well..., independent. Any value other
//! than that is the index of the equation where this dof is the dependent one.
std::vector<int> GetStatusOfDofNumber(const Constraint::Constraints& constraints, DofType dof, int numDofs)
{
    std::vector<int> dofStatus(numDofs, INDEPENDENT); // initialize as INDEPENDENT
    for (int iEquation = 0; iEquation < constraints.GetNumEquations(dof); ++iEquation)
    {
        int dependentDofNumber = constraints.GetEquation(dof, iEquation).GetDependentDofNumber();
        dofStatus.at(dependentDofNumber) = iEquation;
    }
    return dofStatus;
}

DofInfo DofNumbering::Build(const Group<DofNode>& dofNodes, DofType dof, const Constraint::Constraints& constraints)
{
    for (int iEquation = 0; iEquation < constraints.GetNumEquations(dof); ++iEquation)
        for (const auto& term : constraints.GetEquation(dof, iEquation).GetTerms())
            if (not dofNodes.Contains(term.GetNode()))
                throw Exception(__PRETTY_FUNCTION__,
                                "The constraints for dof " + dof.GetName() +
                                        " contain nodes that are not included in _dofNodes_. You "
                                        "may have selected the wrong nodes (coordinate nodes?) "
                                        "from the mesh.");

    int numDependentDofs = constraints.GetNumEquations(dof);
    const int numDofs = InitialUnconstrainedNumbering(dofNodes);
    std::vector<int> dofStatus = GetStatusOfDofNumber(constraints, dof, numDofs);
    // It is _very_ important that the dependent dofs are ordered _exactly_ like the equations. Only this feature will
    // produce an identity matrix for T2.
    // Equation0 --> dependentDofOfEquation0 = numIndependentDofs + 0
    // Equation1 --> dependentDofOfEquation1 = numIndependentDofs + 1
    // Equation2 --> dependentDofOfEquation2 = numIndependentDofs + 2
    //
    // That is why the equation number is stored in the dofStatus vector.

    int countIndependentDofs = 0;
    const int numIndependentDofs = numDofs - numDependentDofs;
    for (auto& node : dofNodes)
        for (int iComponent = 0; iComponent < node.GetNumValues(); ++iComponent)
        {
            int dofNumber = node.GetDofNumber(iComponent);
            if (dofStatus.at(dofNumber) == INDEPENDENT)
            {
                node.SetDofNumber(iComponent, countIndependentDofs++);
            }
            else
            {
                int equationIndex = dofStatus[dofNumber];
                node.SetDofNumber(iComponent, numIndependentDofs + equationIndex);
            }
        }

    DofInfo dofInfo;
    dofInfo.numDependentDofs[dof] = numDependentDofs;
    dofInfo.numIndependentDofs[dof] = numDofs - numDependentDofs;
    return dofInfo;
}

std::vector<int> DofNumbering::Get(const Group<DofNode>& dofNodes, int component)
{
    std::vector<int> dofNumbers;
    dofNumbers.reserve(dofNodes.Size());
    for (const auto& node : dofNodes)
        dofNumbers.push_back(node.GetDofNumber(component));
    return dofNumbers;
}
