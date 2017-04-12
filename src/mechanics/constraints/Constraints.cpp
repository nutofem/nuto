#include <ostream>
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

void Constraint::Constraints::Add(Node::eDof dof, Equation equation)
{
    mEquations[dof].push_back(equation);
    mConstraintsChanged = true;
}

void Constraint::Constraints::Add(Node::eDof dof, std::vector<Equation> equations)
{
    Equations& dofEquations = mEquations[dof];
    dofEquations.insert(dofEquations.begin(), equations.begin(), equations.end());
    mConstraintsChanged = true;
}

void Constraint::Constraints::RemoveAll()
{
    mEquations.clear();
    mConstraintsChanged = true;
}


void Constraint::Constraints::SetHaveChanged(bool value)
{
    mConstraintsChanged = value;
}

bool Constraint::Constraints::HaveChanged() const
{
    return mConstraintsChanged;
}

Eigen::VectorXd Constraint::Constraints::GetRhs(Node::eDof dof, double time) const
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

SparseMatrixCSRVector2General<double> Constraint::Constraints::BuildConstraintMatrix(Node::eDof dof, int nDofs) const
{
    const auto it = mEquations.find(dof);
    if (it == mEquations.end())
        return SparseMatrixCSRVector2General<double>(0, nDofs); // no equations for this dof type

    const Equations& equations = it->second;
    int numEquations = equations.size();

    SparseMatrixCSRVector2General<double> matrix(numEquations, nDofs);

    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
    {
        const auto& equation = equations[iEquation];
        for (const auto& term : equation.GetTerms())
        {
            if (term.GetNode().GetNum(dof) < term.GetComponent())
            {
                // TODO
                // Currently, equations have no knowledge about dofs and cannot call GetNum(dof). This
                // can only be done here in the constraint. 
                // After adding the separation of dofs, each node has only one dof type and the GetNum() 
                // function can be called (e.g.) in the ctor of term.
                std::ostringstream message;
                message << "Cannot access component " << term.GetComponent() << " from node " << term.GetNode() << ".";
                throw MechanicsException(__PRETTY_FUNCTION__, message.str());
            }

            double coefficient = term.GetCoefficient();
            int globalDofNumber = term.GetNode().GetDof(dof, term.GetComponent());
            if (std::abs(coefficient) > 1.e-18)
                matrix.AddValue(iEquation, globalDofNumber, coefficient);
        }
    }
    return matrix;
}

int Constraint::Constraints::GetNumEquations(Node::eDof dof) const
{
    auto it = mEquations.find(dof);
    if (it == mEquations.end())
        return 0;

    return it->second.size();
}

void Constraint::Constraints::ExchangeNodePtr(const NodeBase& oldNode, const NodeBase& newNode)
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
    constexpr double rhsEvaluateTime0 = 0;
    constexpr double rhsEvaluateTime1 = 1;
    for (const auto& dofEquationPair : constraints.mEquations)
    {
        out << "Equations for dof type " << Node::DofToString(dofEquationPair.first) << ":\n";
        const auto& equations = dofEquationPair.second;
        for (const auto& equation : equations)
        {
            out << "Rhs [ " << equation.GetRhs(rhsEvaluateTime0) << " ... " << equation.GetRhs(rhsEvaluateTime1)
                << " ] = ";
            for (const auto& term : equation.GetTerms())
            {
                out << term.GetCoefficient() << " * (dof "
                    << term.GetNode().GetDof(dofEquationPair.first, term.GetComponent()) << " [component "
                    << term.GetComponent() << "]) + ";
            }
            out << "\b\b  \n"; // overrides the last "+" by moving the cursor left twice(\b) and override with blanks
        }
    }
    return out;
}
} /* Constraint */
} /* NuTo */
