#pragma once
#include <tuple>

namespace NuTo
{

//! @brief Evaluates the results of a single newton raphson iteration
//! @Tparam TFunction ... Objects that implements the NuTo::Residual interface
template <typename TFunction>
class LineSearch
{
    using ResidualType = typename TFunction::ResidualType;
    using NormType = typename TFunction::NormType;

public:
    LineSearch(NormType tolerance)
        : mTolerance(tolerance)
    {
    }

    //! @brief evaluates the inputs based on the TUseLinesearch option
    //! @tparam TUseLinesearch ... to use or not to use
    //! @param x ... previous value of the argument x
    //! @param dx ... delta x from the solve
    //! @param f ... function that hat R(x) and Norm(R)
    //! @return ... tuple containing [tolerance reached?, new residual, new x]. These various return values are
    //! an optimization to avoid duplicate calculations
    template <bool TUseLinesearch>
    std::tuple<bool, ResidualType, ResidualType> Evaluate(ResidualType x, ResidualType dx, TFunction f) const
    {
        if (TUseLinesearch)
            return With(x, dx, f);
        else
            return Without(x, dx, f);
    }

private:

    //! @brief calculates the new residual R(x - dx) and compares its norm to the tolerance
    //! @params see Evaluate()
    std::tuple<bool, ResidualType, ResidualType> Without(ResidualType x, ResidualType dx, TFunction f) const
    {
        x -= dx;
        ResidualType r = f.R(x);
        return std::make_tuple(f.Norm(r) < mTolerance, r, x);
    }

    //! @param actually performs linesearch. Applies the x - alpha*dx with decreasing alpha in (0,1)
    //! to ensure a quadratic convergence of newton iteration
    //! @params see Evaluate()
    std::tuple<bool, ResidualType, ResidualType> With(ResidualType x, ResidualType dx, TFunction f) const
    {
        constexpr double minLineSearchStep = 0.01; // 0.5**6 > 0.01 > 0.5**7. Stops after 6 line search steps
        double alpha = 1.;
        NormType norm = f.Norm(f.R(x));
        while (true)
        {
            ResidualType xTrial = x - alpha * dx;
            ResidualType rTrial = f.R(xTrial);
            NormType normTrial = f.Norm(rTrial);

            if (normTrial < mTolerance)
                return std::make_tuple(true, rTrial, xTrial);

            alpha *= 0.5;
            f.Info(alpha, xTrial, rTrial);
            if (alpha < minLineSearchStep)
                return std::make_tuple(false, rTrial, xTrial);
            if (normTrial < (1. - alpha) * norm)
                return std::make_tuple(false, rTrial, xTrial);
        }
    }

    NormType mTolerance;
};
} /* NuTo */
