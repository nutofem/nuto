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

    virtual void Save(const std::string& filename, std::string rType) const
    {
        throw Exception("NuTo::SparseDirectSolver::Save] To be implemented.");
    }

    virtual void Restore(const std::string& filename, std::string rType)
    {
        throw Exception("NuTo::SparseDirectSolver::Restore] To be implemented.");
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
