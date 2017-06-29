#include "math/CubicSplineInterpolation.h"
#include "base/Exception.h"

NuTo::Math::CubicSplineInterpolation::CubicSplineInterpolation(std::vector<std::array<double, 2>> data) : Interpolation::Interpolation(data, 2)
{
    const unsigned n = mData.size();
    ddy.reserve(n);
    std::vector<double> u;
    u.reserve(n - 1);
    // natural boundary condition
    ddy[0] = u[0] = 0.0;
    ddy[n-1] = 0.0;
    for (unsigned i = 1; i < n - 1; ++i)
    {
        double sigma = (mData[i][0] - mData[i-1][0]) / (mData[i+1][0] - mData[i-1][0]);
        double p = sigma*ddy[i-1] + 2.0;
        ddy[i] = (sigma - 1.0) / p;
        u[i]=(mData[i+1][1] - mData[i][1]) / (mData[i+1][0] - mData[i][0])
           - (mData[i][1] - mData[i-1][1]) / (mData[i][0] - mData[i-1][0]);
        u[i]=(6.0*u[i]/(mData[i+1][0] - mData[i-1][0]) - sigma*u[i-1])/p;
    }
    // backsubstitution
    for (unsigned k = n - 2; k > 0; --k)
    {
        ddy[k] = ddy[k]*ddy[k+1] + u[k];
    }
}

double NuTo::Math::CubicSplineInterpolation::operator()(double x)
{
    if (x < mData[0][0] or x > mData.back()[0])
        throw NuTo::Exception("Input x is not within data range of supplied array");

    unsigned index = bisection(x);
    double h = mData[index+1][0] - mData[index][0];
    double a = (mData[index+1][0] - x) / h;
    double b = (x - mData[index][0]) / h;
    double y = a*mData[index][1] + b*mData[index+1][1]
             + ((a*a*a - a) * ddy[index] + (b*b*b - b) * ddy[index+1]) * (h*h) / 6.0;
    return y;
}

double NuTo::Math::CubicSplineInterpolation::derivative(double x)
{
    if (x < mData[0][0] or x > mData.back()[0])
        throw NuTo::Exception("Input x is not within data range of supplied array");

    unsigned index = bisection(x);
    double m = (mData[index+1][1] - mData[index][1]) / (mData[index+1][0] - mData[index][0]);
    return m;
}
