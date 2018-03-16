#include "nuto/base/Timer.h"

#include <iostream>
#include <algorithm>
#include <iomanip>

#include "nuto/base/Logger.h"

#ifdef _OPENMP
#include <ctime>
#endif // _OPENMP


NuTo::Timer::Timer(std::string rMsg, bool rShowTime)
    : mMsg(rMsg)
    , mShowTime(rShowTime)
    , mLogger(nullptr)
#ifdef _OPENMP
    , mCPUTimeInit(clock())
#endif // _OPENMP
    , mWallTimeInit(std::chrono::system_clock::now())
{
}

NuTo::Timer::Timer(std::string rMsg, bool rShowTime, Logger& rLogger)
    : mMsg(rMsg)
    , mShowTime(rShowTime)
    , mLogger(&rLogger)
#ifdef _OPENMP
    , mCPUTimeInit(clock())
#endif // _OPENMP
    , mWallTimeInit(std::chrono::system_clock::now())
{
}

NuTo::Timer::~Timer()
{
    Reset();
}

void NuTo::Timer::Reset()
{
    if (mShowTime)
    {
#ifdef _OPENMP
        double cpuTimeDifference = GetCPUTimeDifference();
#endif // _OPENMP
        double wallTimeDifference = GetTimeDifference();

        std::ostringstream timing;
        timing << "W:" << std::scientific << std::setprecision(2) << wallTimeDifference << "s";
#ifdef _OPENMP
        timing << "  C:" << cpuTimeDifference << "s";
        timing << "  S:" << std::fixed << cpuTimeDifference / wallTimeDifference;
#endif // _OPENMP


        int lengthMsg = static_cast<int>(mMsg.length());
        int lengthTime = static_cast<int>(timing.str().length());


        int numAdditionalBlanks = std::max(0, mMinOutputLength - lengthMsg - lengthTime - 2); // -2 for []
        std::ostringstream out;
        out << "[" << mMsg << "]" << std::string(numAdditionalBlanks, '.') << timing.str() << '\n';
        if (mLogger == nullptr)
            std::cout << out.str();
        else
            *mLogger << out.str();
    }
    mWallTimeInit = std::chrono::system_clock::now();
#ifdef _OPENMP
    mCPUTimeInit = clock();
#endif // _OPENMP
}

void NuTo::Timer::Reset(std::string rMsg)
{
    Reset();
    mMsg = rMsg;
}

//! @brief returns the time from ctor to now in seconds
double NuTo::Timer::GetTimeDifference() const
{
    const auto timeDiff = std::chrono::system_clock::now() - mWallTimeInit;
    return std::chrono::duration_cast<std::chrono::duration<double>>(timeDiff).count();
}

#ifdef _OPENMP
//! @brief returns the time from ctor to now in seconds
double NuTo::Timer::GetCPUTimeDifference() const
{
    return difftime(clock(), mCPUTimeInit) / CLOCKS_PER_SEC;
}
#endif // _OPENMP
