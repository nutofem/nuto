#include <algorithm>
#include "nuto/math/Interpolation.h"
#include "nuto/math/MathException.h"

bool CompareDataPairs(const std::array<double, 2>& x1, const std::array<double, 2>& x2)
{
    return x1[0] < x2[0] ? true : false;
}

NuTo::Math::Interpolation::Interpolation(std::vector<std::array<double, 2>> data,
        unsigned interpolationOrder) : mData{data}, mInterpolationOrder{interpolationOrder}
{
    if (mData.size() < interpolationOrder+1)
    {
        throw NuTo::invalid_argument("Input array does not have enough entries to interpolate.");
    }
    std::sort(mData.begin(), mData.end(), CompareDataPairs);
}

unsigned NuTo::Math::Interpolation::bisection(double x)
{
    unsigned lower = 0;
    unsigned upper = mData.size() - 1;
    unsigned pivot;
    while (upper - lower > 1)
    {
        pivot = (upper + lower) / 2; 
        if (x > mData[pivot][0])
        {
            lower = pivot;
        }
        else
        {
            upper = pivot;
        }
    }
    return lower - (mInterpolationOrder - 1)/2;
}
