#pragma once

#include "base/Exception.h"

namespace NuTo
{
class SparseDirectSolver
{
public:
    //! @brief ... default constructor
    SparseDirectSolver()
        : mShowTime(true)
    {
    }

    bool GetShowTime() const
    {
        return mShowTime;
    }

    void SetShowTime(bool showTime)
    {
        mShowTime = showTime;
    }

    virtual ~SparseDirectSolver() = default;

private:
    bool mShowTime;
};
} // namespace NuTo
