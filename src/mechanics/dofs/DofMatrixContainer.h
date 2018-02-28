#pragma once

#include <map>
#include <vector>
#include <eigen3/Eigen/Core>
#include "mechanics/dofs/DofType.h"

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
template <typename T>
class DofMatrixContainer
{
public:
    T& operator()(DofType d0, DofType d1)
    {
        return this->mData[std::make_pair(d0, d1)];
    }

    const T& operator()(DofType d0, DofType d1) const
    {
        return this->mData.at(std::make_pair(d0, d1));
    }

    //! @brief performs _uninitialized addition_ that resizes the data to the length of `rhs`
    DofMatrixContainer& operator+=(const DofMatrixContainer& rhs)
    {
        for (auto& entry : rhs.mData)
        {
            if (mData.find(entry.first) != mData.end())
                mData[entry.first] += entry.second;
            else
                mData[entry.first] = entry.second;
        }
        return *this;
    }

    DofMatrixContainer& operator*=(double scalar)
    {
        for (auto& entry : mData)
            entry.second *= scalar;
        return *this;
    }

    friend DofMatrixContainer operator+(DofMatrixContainer lhs, const DofMatrixContainer& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    friend DofMatrixContainer operator*(DofMatrixContainer lhs, double scalar)
    {
        lhs *= scalar;
        return lhs;
    }

    //! @brief calculates (*this) += `rhs` * `scalar`
    //! @remark This is the common case in the numerical integration - summing up all integration point contributions
    //! scaled with DetJ and the integration point weight.
    void AddScaled(const DofMatrixContainer& rhs, double scalar)
    {
        for (auto& entry : rhs.mData)
        {
            if (mData.find(entry.first) != mData.end())
                mData[entry.first] += entry.second * scalar;
            else
                mData[entry.first] = entry.second * scalar;
        }
    }

    friend std::ostream& operator<<(std::ostream& out, const DofMatrixContainer<T>& dofMatrix)
    {
        for (auto& entry : dofMatrix.mData)
        {
            auto& dofs = entry.first;
            out << "=== " << dofs.first.GetName() << " " << dofs.second.GetName() << " ===" << std::endl;
            out << entry.second << std::endl;
        }
        out << "====" << std::endl;
        return out;
    }

    std::vector<DofType> DofTypes() const
    {
        std::vector<DofType> dofTypes;
        for (const auto& data : mData)
        {
            AddUnique(&dofTypes, data.first.first);
            AddUnique(&dofTypes, data.first.second);
        }
        return dofTypes;
    }

private:
    static void AddUnique(std::vector<DofType>* dofTypes, DofType dofType)
    {
        bool isNew = true;
        for (const auto& d : *dofTypes)
            if (d.Id() == dofType.Id())
            {
                isNew = false;
                break;
            }
        if (isNew)
            dofTypes->push_back(dofType);
    }

    using DofPair = std::pair<DofType, DofType>;
    struct CompareDofPairs
    {
        bool operator()(const DofPair& a, const DofPair& b) const
        {
            // check for equiv(a.first, b.first)
            if (a.first.Id() == b.first.Id())
                return a.second.Id() < b.second.Id();
            else
                return a.first.Id() < b.first.Id();
        }
    };
    std::map<DofPair, T, CompareDofPairs> mData;
};
} /* NuTo */
