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
        return Eigen::SparseVector<float>(numDofs); // no equations for this dof type

    const Equations& equations = mEquations[dof];
    int numEquations = equations.size();
    Eigen::SparseVector<float> globalVector(numDofs);
    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
    {
        globalVector.coeffRef(equations[iEquation].GetDependentDofNumber()) = equations[iEquation].GetRhs(time);
    }
    return globalVector;
}

JKNumbering Constraints::GetJKNumbering(DofType dof, int numDofs) const
{
    int numJ = numDofs - GetNumEquations(dof);
    int numK = GetNumEquations(dof);
    Eigen::VectorXi independentGlobalNumbering(numJ);
    Eigen::VectorXi dependentGlobalNumbering(numK);

    std::vector<bool> isDofConstrained(numDofs, false);
    for (int i = 0; i < numK; i++)
    {
        int globalDofNumber = GetEquation(dof, i).GetDependentDofNumber();
        if (globalDofNumber == -1 /* should be NodeSimple::NOT_SET */)
            throw Exception(__PRETTY_FUNCTION__,
                            "There is no dof numbering for a node in equation" + std::to_string(i) + ".");
        if (globalDofNumber >= numDofs)
            throw Exception(__PRETTY_FUNCTION__,
                            "The provided dof number of the dependent term exceeds "
                            "the total number of dofs in equation " +
                                    std::to_string(i) + ".");
        isDofConstrained[globalDofNumber] = true;
        dependentGlobalNumbering(i) = globalDofNumber;
    }

    int independentCount = 0;
    for (int i = 0; i < numDofs; i++)
    {
        if (isDofConstrained[i])
            continue;
        independentGlobalNumbering(independentCount) = i;
        independentCount++;
    }
    Eigen::VectorXi joinedNumbering(numDofs);
    joinedNumbering << independentGlobalNumbering, dependentGlobalNumbering;

    return JKNumbering(joinedNumbering, numK);
}

Eigen::SparseMatrix<double> Constraints::BuildUnitConstraintMatrix(DofType dof, int numDofs) const
{
    if (not mEquations.Has(dof))
    {
        Eigen::SparseMatrix<double> unitMatrix(numDofs, numDofs);
        unitMatrix.setIdentity();
        return unitMatrix; // no equations for this dof type
    }

    Eigen::VectorXi jkNumbering = GetJKNumbering(dof, numDofs).mIndices;
    Eigen::VectorXi reverseJKNumbering =
            ((Eigen::PermutationMatrix<Eigen::Dynamic>(jkNumbering)).transpose()).eval().indices();

    const Equations& equations = mEquations[dof];
    int numEquations = equations.size();
    int numIndependentDofs = numDofs - numEquations;

    Eigen::SparseMatrix<double> matrix(numDofs, numIndependentDofs);
    matrix.setZero();

    // add unit entry for all independent dofs
    std::vector<Eigen::Triplet<double>> tripletList;
    for (auto i = 0; i < numIndependentDofs; i++)
    {
        tripletList.push_back(Eigen::Triplet<double>(jkNumbering(i), i, 1.));
    }

    // add entries for constraint dofs
    for (int iEquation = 0; iEquation < numEquations; ++iEquation)
    {
        const auto& equation = equations[iEquation];
        for (Term term : equation.GetIndependentTerms())
        {
            double coefficient = term.GetCoefficient();
            int globalDofNumber = term.GetConstrainedDofNumber();
            int independentDofNumber = reverseJKNumbering[globalDofNumber];
            assert(independentDofNumber < numIndependentDofs);
            if (globalDofNumber == -1 /* should be NodeSimple::NOT_SET */)
                throw Exception(__PRETTY_FUNCTION__,
                                "There is no dof numbering for a node in equation" + std::to_string(iEquation) + ".");

            if (std::abs(coefficient) > 1.e-18)
                tripletList.push_back(
                        Eigen::Triplet<double>(equation.GetDependentDofNumber(), independentDofNumber, -coefficient));
        }
    }
    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return matrix;
}

Eigen::SparseMatrix<double> Constraints::BuildUnitConstraintMatrixInv(DofType dof, int numDofs) const
{
    if (not mEquations.Has(dof))
    {
        Eigen::SparseMatrix<double> unitMatrix(numDofs, numDofs);
        unitMatrix.setIdentity();
        return unitMatrix; // no equations for this dof type
    }

    Eigen::VectorXi jkNumbering = GetJKNumbering(dof, numDofs).mIndices;

    const Equations& equations = mEquations[dof];
    int numEquations = equations.size();
    int numIndependentDofs = numDofs - numEquations;

    Eigen::SparseMatrix<double> matrix(numIndependentDofs, numDofs);
    matrix.setZero();

    // add unit entry for all independent dofs
    std::vector<Eigen::Triplet<double>> tripletList;
    for (auto i = 0; i < numIndependentDofs; i++)
    {
        tripletList.push_back(Eigen::Triplet<double>(i, jkNumbering(i), 1.));
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
    Term newDependentTerm = e.GetDependentTerm();

    if (Contains(mDependentTerms, newDependentTerm))
        throw Exception(__PRETTY_FUNCTION__, "The dependent dof of the new equation "
                                             "is already constrained as a dependent dof in another equation.");

    if (Contains(mIndependentTerms, newDependentTerm))
        throw Exception(__PRETTY_FUNCTION__, "The dependent dof of the new equation "
                                             "is already constrained in another equation.");


    // Any term in the new equation must not constrain the same dof as any existing _dependent_ terms
    for (auto& newTerm : e.GetIndependentTerms())
    {
        if (Contains(mDependentTerms, newTerm))
            throw Exception(__PRETTY_FUNCTION__, "One of the independent terms is already constrained "
                                                 "as a dependent dof in another equation");
    }

    mDependentTerms.insert(e.GetDependentTerm());
    for (auto& newTerm : e.GetIndependentTerms())
        mIndependentTerms.insert(newTerm);
}
