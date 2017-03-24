#pragma once
#include "mechanics/nodes/NodeBase.h"
#include "math/SparseMatrix.h"

namespace NuTo
{
namespace Constaint
{

class Term
{
public:
    Term(const NodeBase& node, int component, double coefficient)
        : mNode(node)
        , mComponent(component), mCoefficient(coefficient)
    {
    }

    const NodeBase& GetNode() const
    {
        return mNode;
    }

    int GetComponent() const
    {
        return mComponent;
    }

    double GetCoefficient() const
    {
        return mCoefficient;
    }

private:
    const NodeBase& mNode;
    int mComponent;
    double mCoefficient;
};

class Equation
{
public:
    Equation(double rhs)
        : mRhs(rhs)
    {
    }

    void AddTerm(Term term)
    {
        mTerms.push_back(term);
    }

    const std::vector<Term>& GetTerms() const
    {
        return mTerms;
    }

    double GetRhs() const
    {
        return mRhs;
    }

private:
    double mRhs;
    std::vector<Term> mTerms;
};

class Constraints
{
    using Equations = std::vector<Equation>;

public:
    void AddEquation(NuTo::Node::eDof dof, Equation equation)
    {
        mEquations[dof].push_back(equation);
    }

    Eigen::VectorXd GetRhs(NuTo::Node::eDof dof) const
    {
        const Equations& equations = mEquations.at(dof);
        int numEquations = equations.size();
        Eigen::VectorXd rhs(numEquations);
        for (int iEquation = 0; iEquation < numEquations; ++iEquation)
            rhs[iEquation] = equations[iEquation].GetRhs();
        return rhs;
    }

    void GetConstraintMatrix(NuTo::SparseMatrix<double>& rConstraintMatrix, NuTo::Node::eDof dof) const
    {
        const Equations& equations = mEquations.at(dof);
        int numEquations = equations.size();

        for (int iEquation = 0; iEquation < numEquations; ++iEquation)
        {
            const auto& equation = equations[iEquation];
            for (const auto& term : equation.GetTerms())
            {
                double coefficient = term.GetCoefficient();
                int globalDofNumber = term.GetNode().GetDof(dof, term.GetComponent());
                if (std::abs(coefficient > 1.e-18))
                    rConstraintMatrix.AddValue(iEquation, globalDofNumber, coefficient);
            }
        }
    }

private:
    std::map<NuTo::Node::eDof, Equations> mEquations;
};

} /* Constaint */
} /* NuTo */
