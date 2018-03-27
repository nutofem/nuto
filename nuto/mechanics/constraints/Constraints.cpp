#include "nuto/mechanics/constraints/Constraints.h"

using namespace NuTo;
using namespace NuTo::Constraint;


void Constraints::Add(DofType dof, Equation equation)
{
    mTermChecker.CheckEquation(equation);

    mEquations[dof].push_back(equation);
}

void Constraints::Add(DofType dof, std::vector<Equation> equations)
{
    if (equations.empty())
        throw Exception(__PRETTY_FUNCTION__, "You called this method with an empty collection of equations. Maybe "
                                             "check the node selection for empty node groups.");
    for (auto equation : equations)
        Add(dof, equation);
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

Eigen::SparseVector<double> Constraints::GetSparseGlobalRhs(DofType dof, int numDofs, double time) const
{
    if (not mEquations.Has(dof))
        return Eigen::SparseVector<float>(0); // no equations for this dof type

    const Equations& equations = mEquations[dof];
    int numEquations = equations.size();
    Eigen::SparseVector<float> globalVector(numDofs);
    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
    {
        globalVector.coeffRef(equations[iEquation].GetDependentDofNumber()) = equations[iEquation].GetRhs(time);
    }
    return globalVector;
}

Eigen::SparseMatrix<double> Constraints::BuildUnitConstraintMatrix2(DofType dof, int numDofs) const
{
    if (not mEquations.Has(dof))
    {
        Eigen::SparseMatrix<double> unitMatrix(numDofs,numDofs);
        unitMatrix.setIdentity();
        return unitMatrix; // no equations for this dof type
    }

    const Equations& equations = mEquations[dof];
    int numEquations = equations.size();
    int numIndependentDofs = numDofs - numEquations;

    Eigen::SparseMatrix<double> matrix(numDofs, numIndependentDofs);
    matrix.setZero();
    std::vector<Eigen::Triplet<double> > tripletList;
    for (auto i=0; i<numIndependentDofs; i++)
    {
        tripletList.push_back(Eigen::Triplet<double>(i, i, 1.));
    }

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
                                                         "is not in accordance with the equation numbering.");
                continue; // Do not put the value into the matrix.
            }

            if (std::abs(coefficient) > 1.e-18)
                tripletList.push_back(Eigen::Triplet<double>(numIndependentDofs + iEquation, globalDofNumber, -coefficient));
        }
    }
    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
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

template <typename T>
bool Contains(const T& container, NuTo::Constraint::Term t)
{
    return container.find(t) != container.end();
}

bool Constraints::TermChecker::TermCompare::operator()(const Term& lhs, const Term& rhs) const
{
    if (&lhs.GetNode() != &rhs.GetNode())
        return &lhs.GetNode() < &rhs.GetNode();
    return lhs.GetComponent() < rhs.GetComponent();
}

void Constraints::TermChecker::CheckEquation(Equation e)
{
    // The new dependent term must not constrain the same dof as any existing terms
    Term newDependentTerm = e.GetTerms()[0];

    if (Contains(mDependentTerms, newDependentTerm))
        throw Exception(__PRETTY_FUNCTION__, "The dependent dof of the new equation "
                                             "is already constrained as a dependent dof in another equation.");

    if (Contains(mOtherTerms, newDependentTerm))
        throw Exception(__PRETTY_FUNCTION__, "The dependent dof of the new equation "
                                             "is already constrained in another equation.");


    // Any term in the new equation must not constrain the same dof as any existing _dependent_ terms
    // since the dependent term of the new equation is already checked above, we omit it here and start loop at 1.
    for (size_t iNewTerm = 1; iNewTerm < e.GetTerms().size(); ++iNewTerm)
    {
        Term newTerm = e.GetTerms()[iNewTerm];
        if (Contains(mDependentTerms, newTerm))
            throw Exception(__PRETTY_FUNCTION__, "One of the new terms is already constrained "
                                                 "as a dependent dof in another equation");
    }

    mDependentTerms.insert(e.GetTerms()[0]);
    for (size_t iNewTerm = 1; iNewTerm < e.GetTerms().size(); ++iNewTerm)
        mOtherTerms.insert(e.GetTerms()[iNewTerm]);
}
