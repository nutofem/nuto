#pragma once

#include "nuto/math/Equation.h"
#include "ConstraintApplicator.h"

namespace NuTo
{

class LagrangeApplicator : public ConstraintApplicator
{
public:
    LagrangeApplicator(std::vector<Math::Equation> eqs)
        : mEquations(eqs)
    {
    }

    Eigen::SparseMatrix<double> Apply(Eigen::SparseMatrix<double>& K) override
    {
        size_t originalRows = K.rows();
        size_t newSize = originalRows + mEquations.size();
        K.conservativeResize(newSize, newSize);
        size_t i = 0;
        for (auto eq : mEquations)
        {
            for (auto term : eq.Terms())
            {
                K.insert(originalRows + i, term.Id()) = term.Coefficient();
                K.insert(term.Id(), originalRows + i) = term.Coefficient();
            }
            i++;
        }
        K.makeCompressed();

        return K;
    }

    Eigen::VectorXd Apply(Eigen::VectorXd& f) override
    {
        size_t originalRows = f.rows();
        size_t newSize = originalRows + mEquations.size();
        f.conservativeResize(newSize);
        size_t i = 0;
        for (auto eq : mEquations)
        {
            f[originalRows + i] = eq.Value();
            i++;
        }
        return f;
    }

    Eigen::VectorXd ConvertBack(Eigen::VectorXd& u) override
    {
        size_t newSize = u.rows() - mEquations.size();
        return u.head(newSize);
    }

private:
    std::vector<Math::Equation> mEquations;
};

} // namespace NuTo
