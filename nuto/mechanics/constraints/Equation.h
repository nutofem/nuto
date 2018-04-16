#pragma once
#include <vector>
#include "nuto/mechanics/constraints/Term.h"

namespace NuTo
{
namespace Constraint
{
typedef std::function<double(double)> RhsFunction;

//! @brief stores a constraint equation
class Equation
{
public:
    //! @brief ctor with constant rhs, defaults to 0
    //! @param dependentNode node reference
    //! @param dependentComponent component in the dof vector of the node
    //! @param rhs value for the constant rhs
    Equation(const DofNode& dependentNode, int dependentComponent, RhsFunction rhs)
        : mRhs(rhs)
        , mDependentTerm(dependentNode, dependentComponent, 1)
    {
    }

    //! @brief adds a term to the equation
    //! @param term equation term
    void AddIndependentTerm(Term term)
    {
        mIndependentTerms.push_back(term);
    }

    //! @brief evaluates the rhs function at a given time
    //! @param time global time
    //! @return rhs(time)
    double GetRhs(double time) const
    {
        return mRhs(time);
    }

    //! @brief getter for mTerms
    const Term& GetDependentTerm() const
    {
        return mDependentTerm;
    }

    //! @brief getter for mTerms
    const std::vector<Term>& GetIndependentTerms() const
    {
        return mIndependentTerms;
    }

    int GetDependentDofNumber() const
    {
        return mDependentTerm.GetConstrainedDofNumber();
    }

private:
    //! @brief rhs function
    RhsFunction mRhs;

    //! @brief dependent term rhs(t) = dependentTerm + independentTerms
    Term mDependentTerm;

    //! @brief all other terms
    std::vector<Term> mIndependentTerms;
};

} /* Constaint */
} /* NuTo */
