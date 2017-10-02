#pragma once

#include <array>
#include <vector>
#include "math/Interpolation.h"

namespace NuTo
{
namespace Math
{
class LinearInterpolation : public Interpolation
{
public:
    //! create interpolation object; call with data array
    LinearInterpolation(std::vector<std::array<double, 2>> data)
        : Interpolation{data, 2} {};

    //! return interpolated value at x
    double operator()(double x) const override;

    //! calculate first derivative at x
    double derivative(double x) const override;
};
} // namespace Math
} // namespace NuTo
