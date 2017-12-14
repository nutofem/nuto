#pragma once
#include <vector>
#include "mechanics/constraints/Term.h"

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
    Equation(const NodeSimple& dependentNode, int dependentComponent, RhsFunction rhs)
        : mRhs(rhs)
        , mTerms{Term(dependentNode, dependentComponent, 1)}
    {
    }

    //! @brief adds a term to the equation
    //! @param term equation term
    void AddTerm(Term term)
    {
        mTerms.push_back(term);
    }

    //! @brief evaluates the rhs function at a given time
    //! @param time global time
    //! @return rhs(time)
    double GetRhs(double time) const
    {
        return mRhs(time);
    }

    //! @brief getter for mTerms
    const std::vector<Term>& GetTerms() const
    {
        return mTerms;
    }

    int GetDependentDofNumber() const
    {
        return mTerms.front().GetConstrainedDofNumber();
    }

private:
    //! @brief rhs function
    RhsFunction mRhs;

    //! @brief terms
    std::vector<Term> mTerms;
};

} /* Constaint */
} /* NuTo */
