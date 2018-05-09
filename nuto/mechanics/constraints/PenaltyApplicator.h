#pragma once

#include "nuto/math/Equation.h"
#include "nuto/base/Exception.h"
#include "ConstraintApplicator.h"

namespace NuTo
{

class PenaltyApplicator : public ConstraintApplicator
{
public:
    PenaltyApplicator(std::vector<Math::Equation> eqs)
        : mEquations(eqs)
    {
    }

    Eigen::SparseMatrix<double> Apply(Eigen::SparseMatrix<double>& K) override
    {
        for (auto it = K.valuePtr(); it != K.valuePtr() + K.nonZeros(); ++it)
            if (*it > mPenaltyParameter)
                mPenaltyParameter = *it;
        mPenaltyParameter *= 1e15;

        for (auto eq : mEquations)
        {
            for (auto termI : eq.Terms())
            {
                for (auto termJ : eq.Terms())
                    K.insert(termI.Id(), termJ.Id()) = mPenaltyParameter * termI.Coefficient() * termJ.Coefficient();
            }
        }
        K.makeCompressed();

        return K;
    }

    Eigen::VectorXd Apply(Eigen::VectorXd& f) override
    {
        if (mPenaltyParameter == -1.0)
            throw Exception(__PRETTY_FUNCTION__,
                            "You need to call Apply on a matrix first to calculate the size of the penalty parameter.");

        for (auto eq : mEquations)
        {
            for (auto termI : eq.Terms())
                f[termI.Id()] += mPenaltyParameter * termI.Coefficient() * eq.Value();
        }
        return f;
    }

    Eigen::VectorXd ConvertBack(Eigen::VectorXd& u) override
    {
        return u;
    }

private:
    std::vector<Math::Equation> mEquations;
    double mPenaltyParameter = -1.0;
};

} // namespace NuTo
