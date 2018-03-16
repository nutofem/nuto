#pragma once

#include <map>
#include "nuto/mechanics/dofs/DofType.h"

namespace NuTo
{
template <typename T>
class DofContainer
{
public:
    virtual ~DofContainer() = default;

    //! @brief nonconst access, similar to map::operator[]()
    //! @param dofType dof type
    //! @return reference to either
    //          an existing value
    //          an newly default constructed value
    //! @remark This requires T to be default constructable.
    T& operator[](DofType dofType)
    {
        return mData[dofType.Id()];
    }

    //! @brief nonconst access, similar to map::at()
    //! @param dofType dof type
    //! @return reference to an existing value, throws if there is no value
    //! @remark This does not default construct a new T and thus does not
    //          require T to be default constructable.
    T& At(DofType dofType)
    {
        return mData.at(dofType.Id());
    }

    //! @brief const access
    //! @param dofType dof type
    //! @return const reference to existing value, throws if there is no value
    const T& operator[](DofType dofType) const
    {
        return mData.at(dofType.Id());
    }

    //! @brief copies a `t` into the container, throws, if there already is an entry at `dofType`
    //! @param dofType dof type
    //! @param t value to insert
    void Insert(DofType dofType, T t)
    {
        auto it = mData.emplace(dofType.Id(), t); // it = pair<iterator, bool>
        if (not it.second)
            throw Exception(__PRETTY_FUNCTION__,
                            "Insert failed. Container already contains an entry for " + dofType.GetName() + ".");
    }

    bool Has(DofType dof) const
    {
        return mData.find(dof.Id()) != mData.end();
    }

protected:
    std::map<int, T> mData;
};
} /* NuTo */
