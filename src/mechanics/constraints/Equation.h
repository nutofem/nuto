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
    //! @param rhs value for the constant rhs
    Equation(double rhs = 0)
        : mRhs([=](double) { return rhs; })
    {
    }

    //! @brief ctor with terms and constant rhs
    //! @param terms equation terms
    //! @param rhs value for the constant rhs
    Equation(std::vector<Term> terms, double rhs = 0)
        : mRhs([=](double) { return rhs; })
        , mTerms(terms)
    {
    }

    //! @brief ctor with a rhs function
    //! @param rhs rhs function
    Equation(RhsFunction rhs)
        : mRhs(rhs)
    {
    }

    //! @brief ctor with a rhs function and terms
    //! @param rhs rhs function
    //! @param terms equation terms
    Equation(RhsFunction rhs, std::vector<Term> terms)
        : mRhs(rhs)
        , mTerms(terms)
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

    //! @brief nonconst getter for mTerms
    //! @remark IMO only used for ExchangeNodePtr which is very questionable...
    std::vector<Term>& GetTerms()
    {
        return mTerms;
    }

private:
    //! @brief rhs function
    RhsFunction mRhs;

    //! @brief terms
    std::vector<Term> mTerms;
};

} /* Constaint */
} /* NuTo */
