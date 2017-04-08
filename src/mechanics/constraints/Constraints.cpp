#include <ostream>
#include "mechanics/constraints/Constraints.h"
#include "mechanics/nodes/NodeEnum.h"

void NuTo::Constraint::Constraints::Add(Node::eDof dof, Equation equation)
{
    mEquations[dof].push_back(equation);
    mHasNewConstraints = true;
}

void NuTo::Constraint::Constraints::Add(Node::eDof dof, std::vector<Equation> equations)
{
    Equations& dofEquations = mEquations[dof];
    dofEquations.insert(dofEquations.begin(), equations.begin(), equations.end());
    mHasNewConstraints = true;
}

void NuTo::Constraint::Constraints::SetHasNewConstraints(bool value)
{
    mHasNewConstraints = value;
}

bool NuTo::Constraint::Constraints::HasNewConstraints() const
{
    return mHasNewConstraints;
}

Eigen::VectorXd NuTo::Constraint::Constraints::GetRhs(Node::eDof dof, double time) const
{
    const auto it = mEquations.find(dof);
    if (it == mEquations.end())
        return Eigen::VectorXd::Zero(0); // no equations for this dof type

    const Equations& equations = it->second;
    int numEquations = equations.size();
    Eigen::VectorXd rhs(numEquations);
    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
        rhs[iEquation] = equations[iEquation].GetRhs(time);
    return rhs;
}

void NuTo::Constraint::Constraints::BuildConstraintMatrix(SparseMatrix<double>& rConstraintMatrix, Node::eDof dof) const
{
    const auto it = mEquations.find(dof);
    if (it == mEquations.end())
        return; // no equations for this dof type

    const Equations& equations = it->second;
    int numEquations = equations.size();

    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
    {
        const auto& equation = equations[iEquation];
        for (const auto& term : equation.GetTerms())
        {
            if (term.GetNode().GetNum(dof) < term.GetComponent())
                throw MechanicsException(__PRETTY_FUNCTION__,
                                         "Cannot access component " + std::to_string(term.GetComponent()) +
                                                 " from NuTo::NodeBase " + term.GetNode().GetNodeTypeStr());

            double coefficient = term.GetCoefficient();
            int globalDofNumber = term.GetNode().GetDof(dof, term.GetComponent());
            if (std::abs(coefficient > 1.e-18))
                rConstraintMatrix.AddValue(iEquation, globalDofNumber, coefficient);
        }
    }
}

int NuTo::Constraint::Constraints::GetNumEquations(Node::eDof dof) const
{
    auto it = mEquations.find(dof);
    if (it == mEquations.end())
        return 0;

    return it->second.size();
}

void NuTo::Constraint::Constraints::ExchangeNodePtr(const NodeBase& oldNode, const NodeBase& newNode)
{
    for (auto& mapPair : mEquations)
    {
        auto& equations = mapPair.second;
        for (auto& equation : equations)
        {
            for (auto& term : equation.GetTerms())
                term.ExchangeNode(oldNode, newNode);
        }
    }
}

namespace NuTo
{
namespace Constraint
{
std::ostream& operator<<(std::ostream& out, const Constraints& constraints)
{
    for (const auto& dofEquationPair : constraints.mEquations)
    {
        out << "Equations for dof type " << Node::DofToString(dofEquationPair.first) << ":\n";
        const auto& equations = dofEquationPair.second;
        for (const auto& equation : equations)
        {
            out << "Rhs(TODO) = ";
            for (const auto& term : equation.GetTerms())
            {
                out << term.GetCoefficient() << " * (dof " 
                    << term.GetNode().GetDof(dofEquationPair.first, term.GetComponent()) 
                    << " [comp " << term.GetComponent() << "]) + ";
            }
            out << '\n';
        }
    }
    return out;
}
} /* Constraint */
} /* NuTo */
