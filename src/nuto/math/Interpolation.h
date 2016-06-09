#pragma once

#include <array>
#include <vector>

namespace NuTo {
namespace Math
{
class Interpolation
{
public:

    //! create interpolation object; call with data array
    Interpolation(std::vector<std::array<double, 2>> data);

    //! return interpolated value at x
    virtual double operator()(double x) = 0;

    //! calculate first derivative at x
    virtual double derivative(double x) = 0;

protected:

    std::vector<std::array<double, 2>> mData;

    int bisection(double x);

};
} // namespace Math
} // namespace NuTo
