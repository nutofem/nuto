#include "mechanics/constraintsPde/Constraints.h"

using namespace NuTo;
using namespace NuTo::ConstraintPde;

void Constraints::Add(DofType dof, Equation equation)
{
    mEquations[dof].push_back(equation);
    mConstraintsChanged = true;
}

void Constraints::Add(DofType dof, std::vector<Equation> equations)
{
    Equations& dofEquations = mEquations[dof];
    dofEquations.insert(dofEquations.begin(), equations.begin(), equations.end());
    mConstraintsChanged = true;
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

Eigen::SparseMatrix<double> Constraints::BuildConstraintMatrix(DofType dof, int nDofs) const
{
    if (not mEquations.Has(dof))
        return Eigen::SparseMatrix<double>(0, nDofs); // no equations for this dof type

    const Equations& equations = mEquations[dof]; 
    int numEquations = equations.size();

    Eigen::SparseMatrix<double> matrix(numEquations, nDofs);

    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
    {
        const auto& equation = equations[iEquation];
        for (const auto& term : equation.GetTerms())
        {
            double coefficient = term.GetCoefficient();
            int globalDofNumber = term.GetConstrainedDofNumber();
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

