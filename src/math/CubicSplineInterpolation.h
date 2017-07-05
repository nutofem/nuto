#pragma once

#include <array>
#include <vector>
#include "math/Interpolation.h"

namespace NuTo
{
namespace Math
{
class CubicSplineInterpolation : public Interpolation
{
public:
    //! create interpolation object; call with data array
    CubicSplineInterpolation(std::vector<std::array<double, 2>> data);

    //! return interpolated value at x
    double operator()(double x) override;

    //! calculate first derivative at x
    double derivative(double x) override;

private:
    std::vector<double> ddy;
};
} // namespace Math
} // namespace NuTo
