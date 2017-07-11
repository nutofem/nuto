#include "math/LinearInterpolation.h"
#include "base/Exception.h"

double NuTo::Math::LinearInterpolation::operator()(double x)
{
    if (x < mData[0][0] or x > mData.back()[0])
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Input x is not within data range of supplied array");

    unsigned index = bisection(x);
    double m = (mData[index + 1][1] - mData[index][1]) / (mData[index + 1][0] - mData[index][0]);
    double n = mData[index][1];
    double x_tmp = x - mData[index][0];
    return m * x_tmp + n;
}

double NuTo::Math::LinearInterpolation::derivative(double x)
{
    if (x < mData[0][0] or x > mData.back()[0])
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Input x is not within data range of supplied array");

    unsigned index = bisection(x);
    double m = (mData[index + 1][1] - mData[index][1]) / (mData[index + 1][0] - mData[index][0]);
    return m;
}
