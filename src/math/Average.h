#pragma once
namespace NuTo
{

//! @brief calculates the cumulative moving average
struct Average
{
    long double mCurrentAverage = 0;
    int mNum = 0;

    void operator()(int x)
    {
        mCurrentAverage = (x + mNum * mCurrentAverage)/(mNum + 1);
        mNum++;
    }
};
}
