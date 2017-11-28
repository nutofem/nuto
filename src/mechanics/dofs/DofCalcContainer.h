#pragma once

#include <vector>
#include <eigen3/Eigen/Core>
#include "mechanics/dofs/DofType.h"

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
//! @tparam T type of the data that supports all the arithmetic operations
template <typename T>
class DofCalcContainer
{
public:
    //! @brief initializes the data container zu the most common size 1
    //! @remark explicit to avoid construction with `double` type that caused ambiguities in some overloaded methods
    explicit DofCalcContainer(int size = 1)
        : mData(size)
    {
    }

    //! @brief const access
    //! @param dofType dof type
    //! @return const reference to existing value, throws if there is no value
    const T& operator[](const DofType& dofType) const
    {
        return mData.at(dofType.Id());
    }

    //! @brief nonconst access
    //! @param dofType dof type
    //! @return reference to either
    //!         an existing value
    //!         an newly default constructed value
    //! @remark This requires T to be default constructable.
    T& operator[](const DofType& dofType)
    {
        return ResizingIdAccess(dofType.Id());
    }

    //! @brief performs _uninitialized addition_ that resizes the data to the length of `rhs`
    DofCalcContainer& operator+=(const DofCalcContainer& rhs)
    {
        for (size_t i = 0; i < rhs.mData.size(); ++i)
        {
            auto& thisData = ResizingIdAccess(i);
            if (thisData.size() == 0)
                thisData = rhs.mData[i];
            else
                thisData += rhs.mData[i];
        }
        return *this;
    }

    DofCalcContainer& operator*=(double scalar)
    {
        for (auto& data : this->mData)
            data *= scalar;
        return *this;
    }

    friend DofCalcContainer operator+(DofCalcContainer lhs, DofCalcContainer&& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    friend DofCalcContainer operator*(DofCalcContainer lhs, double scalar)
    {
        lhs *= scalar;
        return lhs;
    }

    friend std::ostream& operator<<(std::ostream& out, const DofCalcContainer<T>& v)
    {
        for (size_t i = 0; i < v.mData.size(); ++i)
        {
            out << "==== Index " << i << "====\n";
            out << v.mData[i] << '\n';
        }
        return out;
    }

    //! @brief calculates (*this) += `rhs` * `scalar`
    //! @remark This is the common case in the numerical integration - summing up all integration point contributions
    //! scaled with DetJ and the integration point weight.
    void AddScaled(const DofCalcContainer& rhs, double scalar)
    {
        for (int i = 0; i < rhs.mData.size(); ++i)
        {
            auto& thisData = ResizingIdAccess(i);
            if (thisData.size() == 0)
                thisData = rhs.mData[i] * scalar;
            else
                thisData += rhs.mData[i] * scalar;
        }
    }


protected:
    //! @brief data container
    std::vector<T> mData;

    //! @brief access to `i`-th entry of the data. Default constructs data up to `i`, if needed.
    //! @param i data index
    //! @return nonconst reference to `i`-th data entry
    T& ResizingIdAccess(int i)
    {
        if (static_cast<int>(mData.size()) < i + 1)
            mData.resize(i + 1);
        return mData[i];
    }
};
} /* NuTo */
