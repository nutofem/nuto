#include "mechanics/constraintsPde/Constraints.h"

using namespace NuTo;
using namespace NuTo::ConstraintPde;

void Constraints::Add(DofType dof, Equation equation)
{
    // The new dependent term must not constrain the same dof as any existing terms
    Term newDependentTerm = equation.GetTerms()[0];
    for (int iExistingEq = 0; iExistingEq < GetNumEquations(dof); ++iExistingEq)
    {
        for (Term existingTerm : GetEquation(dof, iExistingEq).GetTerms())
        {
            if (&newDependentTerm.GetNode() == &existingTerm.GetNode() &&
                newDependentTerm.GetComponent() == existingTerm.GetComponent())
                throw Exception(__PRETTY_FUNCTION__, "The dependent dof of the new equation "
                                                     "is already constrained by equation " +
                                                             std::to_string(iExistingEq) + ".");
        }
    }


    // Any term in the new equation must not constrain the same dof as any existing _dependent_ terms
    // since the dependent term of the new equation is already checked above, we omit it here and start loop at 1.
    for (int iNewTerm = 1; iNewTerm < equation.GetTerms().size(); ++iNewTerm)
    {
        Term newTerm = equation.GetTerms()[iNewTerm];
        for (int iExistingEq = 0; iExistingEq < GetNumEquations(dof); ++iExistingEq)
        {
            Term existingDependentTerm = GetEquation(dof, iExistingEq).GetTerms()[0];
            if (&newTerm.GetNode() == &existingDependentTerm.GetNode() &&
                newTerm.GetComponent() == existingDependentTerm.GetComponent())
                throw Exception(__PRETTY_FUNCTION__, "One of the new terms is already constrained "
                                                     "as a dependent dof in equation " +
                                                             std::to_string(iExistingEq) + ".");
        }
    }

    mEquations[dof].push_back(equation);
    mConstraintsChanged = true;
}

void Constraints::Add(DofType dof, std::vector<Equation> equations)
{
    for (auto equation : equations)
        Add(dof, equation);
}

void Constraints::SetHaveChanged(bool value)
{
    mConstraintsChanged = value;
}

bool Constraints::HaveChanged() const
{
    return mConstraintsChanged;
}

Eigen::VectorXd Constraints::GetRhs(DofType dof, double time) const
{
    if (not mEquations.Has(dof))
        return Eigen::VectorXd::Zero(0); // no equations for this dof type

    const Equations& equations = mEquations[dof];
    int numEquations = equations.size();
    Eigen::VectorXd rhs(numEquations);
    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
        rhs[iEquation] = equations[iEquation].GetRhs(time);
    return rhs;
}

Eigen::SparseMatrix<double> Constraints::BuildConstraintMatrix(DofType dof, int numIndependentDofs) const
{
    if (not mEquations.Has(dof))
        return Eigen::SparseMatrix<double>(0, numIndependentDofs); // no equations for this dof type

    const Equations& equations = mEquations[dof];
    int numEquations = equations.size();

    Eigen::SparseMatrix<double> matrix(numEquations, numIndependentDofs);

    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
    {
        const auto& equation = equations[iEquation];
        for (Term term : equation.GetTerms())
        {
            double coefficient = term.GetCoefficient();
            int globalDofNumber = term.GetConstrainedDofNumber();
            if (globalDofNumber == -1 /* should be NodeSimple::NOT_SET */)
                throw Exception(__PRETTY_FUNCTION__,
                                "There is no dof numbering for a node in equation" + std::to_string(iEquation) + ".");

            if (globalDofNumber >= numIndependentDofs)
            {
                // This corresponds to the last block of size numDependentDofs x numDependentDofs that should be an
                // identity matrix
                int transformedDof = globalDofNumber - numIndependentDofs;
                if (transformedDof != iEquation)
                    throw Exception(__PRETTY_FUNCTION__, "The numbering of the dependent dofs "
                                                         "is not in accordance to the equation numbering.");
                continue; // Do not put the value into the matrix.
            }

            if (std::abs(coefficient) > 1.e-18)
                matrix.coeffRef(iEquation, globalDofNumber) = coefficient;
        }
    }
    return matrix;
}

int Constraints::GetNumEquations(DofType dof) const
{
    if (not mEquations.Has(dof))
        return 0;
    return mEquations[dof].size();
}

const Equation& Constraints::GetEquation(DofType dof, int equationNumber) const
{
    return mEquations[dof].at(equationNumber);
}
