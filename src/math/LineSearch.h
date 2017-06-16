#pragma once
#include <tuple> // for std::pair

namespace NuTo
{

template <typename TFunction>
class LineSearchFalse
{
    using ResidualType = typename TFunction::ResidualType;
    using NormType = typename TFunction::NormType;

public:
    LineSearchFalse(NormType tolerance)
        : mTolerance(tolerance)
    {
    }

    std::tuple<bool, ResidualType, ResidualType> operator()(ResidualType x, ResidualType dx, TFunction f) const
    {
        ResidualType xNew = x - dx;
        ResidualType r = f.R(xNew);
        return std::make_tuple(f.Norm(r) < mTolerance, r, xNew);
    }

private:
    NormType mTolerance;
};

template <typename TFunction>
class LineSearchTrue
{
    using ResidualType = typename TFunction::ResidualType;
    using NormType = typename TFunction::NormType;

public:
    LineSearchTrue(NormType tolerance, double minLineSearchStep)
        : mTolerance(tolerance)
        , mMinLineSearchStep(minLineSearchStep)
    {
    }

    std::tuple<bool, ResidualType, ResidualType> operator()(ResidualType x, ResidualType dx, TFunction f) const
    {
        double alpha = 1.;
        ResidualType r = f.R(x);
        NormType norm = f.Norm(r);
        while (true)
        {
            ResidualType xTrial = x - alpha * dx;
            ResidualType rTrial = f.R(xTrial);
            NormType normTrial = f.Norm(rTrial);

            if (normTrial < mTolerance)
                return std::make_tuple(true, rTrial, xTrial);

            alpha *= 0.5;
            if (alpha < mMinLineSearchStep)
                return std::make_tuple(false, rTrial, xTrial);
            if (normTrial < (1. - alpha) * norm)
                return std::make_tuple(false, rTrial, xTrial);
        }
    }

private:
    NormType mTolerance;
    double mMinLineSearchStep;
};


} /* NuTo */
