#pragma once
#include <vector>
#include "mechanics/constraints/Term.h"

namespace NuTo
{
namespace Constraint
{

typedef std::function<double(double)> RhsFunction;

class Equation
{
public:
    Equation(double rhs = 0)
        : mRhs([=](double){return rhs;})
    {
    }

    Equation(RhsFunction rhs)
        : mRhs(rhs)
    {
    }

    Equation(RhsFunction rhs, std::vector<Term> terms)
        : mRhs(rhs)
        , mTerms(terms)
    {
    }

#ifndef SWIG
    Equation(const Equation&) = default;
    Equation(Equation&&) = default;

    Equation& operator=(const Equation&) = default;
    Equation& operator=(Equation&&) = default;
#endif /* SWIG */

    void AddTerm(Term term)
    {
        mTerms.push_back(term);
    }

    const std::vector<Term>& GetTerms() const
    {
        return mTerms;
    }

    std::vector<Term>& GetTerms()
    {
        return mTerms;
    }

    double GetRhs(double time) const
    {
        return mRhs(time);
    }

    ~Equation() = default;

private:
    RhsFunction mRhs;
    std::vector<Term> mTerms;
};

} /* Constaint */
} /* NuTo */
