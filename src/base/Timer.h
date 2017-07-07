#pragma once

#include <string>
#include "base/Logger.h"

#include <iosfwd>
#include <chrono>
#include <memory>

namespace NuTo
{

//! @brief prints the lifetime of a Timer object on destruction
class Timer
{
public:
    //! @brief ctor, saves the current time
    //! @param rMsg ... msg to print on destruction
    //! @param rShowTime ... false: no output
    Timer(std::string rMsg, bool rShowTime = true);

    //! @brief ctor, saves the current time
    //! @param rMsg ... msg to print on destruction
    //! @param rShowTime ... false: no output
    Timer(std::string rMsg, bool rShowTime, Logger& rLogger);

    Timer(const Timer&) = delete;
    Timer& operator=(const Timer&) = delete;

    //! @brief dtor, prints the msg and the lifetime
    ~Timer();

    void Reset();

    void Reset(std::string rMsg);


    //! @brief returns the time from ctor to now in seconds
    double GetTimeDifference() const;

#ifdef _OPENMP
    //! @brief returns the time from ctor to now in seconds
    double GetCPUTimeDifference() const;
#endif // _OPENMP

    std::string mMsg;
    bool mShowTime;
    Logger* mLogger;

#ifdef _OPENMP
    double mCPUTimeInit;
#endif // _OPENMP

    static constexpr int mMinOutputLength = 90;

    std::chrono::time_point<std::chrono::system_clock> mWallTimeInit;
};
}
