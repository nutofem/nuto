#include <algorithm>
#include "math/Interpolation.h"
#include "math/MathException.h"

bool CompareDataPairs(const std::array<double, 2>& x1, const std::array<double, 2>& x2)
{
    return x1[0] < x2[0] ? true : false;
}

NuTo::Math::Interpolation::Interpolation(std::vector<std::array<double, 2>> data, unsigned numNeighborPoints)
    : mData{data}
    , mNumNeighborPoints{numNeighborPoints}
{
    if (mData.size() < mNumNeighborPoints)
    {
        throw NuTo::invalid_argument("Input array does not have enough entries to interpolate.");
    }
    std::sort(mData.begin(), mData.end(), CompareDataPairs);
}

unsigned NuTo::Math::Interpolation::bisection(double x)
{
    unsigned lower = 0;
    unsigned upper = mData.size() - 1;
    while (upper - lower > 1)
    {
        unsigned pivot = (upper + lower) / 2;
        if (x > mData[pivot][0])
        {
            lower = pivot;
        }
        else
        {
            upper = pivot;
        }
    }
    return lower + (mNumNeighborPoints - 2) / 2;
}
