#pragma once

namespace NuTo
{
//! @brief This class provides a unique id (beginning at 0 and incrementing for each object).
//! @tparam T Without the template, the id would be unique for all classes. Say, 1000 nodes are created first and get
//! the ids 0...999. Element ids would then start at 1000. With this template, T=Node or T=Element, each T starts at 0
//! again.
template <typename T>
class UniqueId
{
public:
    UniqueId()
        : mId(mIdCounter++)
    {
    }

    int Id() const
    {
        return mId;
    }

private:
    int mId;
    static int mIdCounter;
};

template <typename T>
int UniqueId<T>::mIdCounter = 0;

} /* NuTo */
