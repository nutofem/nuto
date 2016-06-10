#pragma once

#include <array>
#include <vector>
#include <functional>

namespace NuTo {
namespace Math
{
class Interpolation
{
public:

    //! create interpolation object; call with data array
    Interpolation(std::vector<std::array<double, 2>> data, unsigned interpolationOrder);

    //! return interpolated value at x; order is the order of the derivative
    virtual double operator()(double x) = 0;

    //! calculate first derivative at x
    virtual double derivative(double x) = 0;

    //! function object
    std::function<double(double)> f = [this](double x) { return operator()(x); };

    //! derivative function object
    std::function<double(double)> df = [this](double x) { return derivative(x); };

protected:

    std::vector<std::array<double, 2>> mData;

    unsigned mInterpolationOrder;

    unsigned bisection(double x);
};
} // namespace Math
} // namespace NuTo
