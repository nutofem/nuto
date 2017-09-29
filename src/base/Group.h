#pragma once

#include <vector>
#include <boost/iterator/indirect_iterator.hpp>

namespace NuTo
{
namespace Groups
{

//! @brief Ordered container class for elements, nodes and the like. No duplicates.
//!
//! The group saves pointers to its members but can be iterated through
//! by value. The members are sorted in order of their insertion.
//!
//! @tparam T type of group members
//! @tparam TCompare functor representing the comparison operation used for sorting
template <typename T>
class Group 
{
    //! @brief Container, similar to std::set, with an underlying std::vector.
    class SortedUniqueVector
    {
    public:
        //! @brief adds an element if it is not already in the container
        //! @param element element to add
        //! @return true if element was added
        bool Add(const T& element)
        {
            auto it = std::lower_bound(mData.begin(), mData.end(), &element);
            if (it == mData.end() || *it != &element)
            {
                mData.insert(it, &element);
                return true;
            }
            return false;
        }

        bool Contains(const T& element) const
        {
            auto it = std::lower_bound(mData.begin(), mData.end(), &element);
            if (it == mData.end())
                return false;
            return *it == &element;
        }

    private:
        std::vector<const T*> mData;
    };

    typedef std::vector<T*> Data;

public:
    //! @brief indirect (dereferencing) iterator to provide value semantics for the iterators
    typedef boost::indirect_iterator<typename Data::iterator> GroupIterator;
    typedef boost::indirect_iterator<typename Data::const_iterator> ConstGroupIterator;

    //! @brief Create an empty group
    Group() = default;

    //! @brief Create a group containing a single element
    //! @param element single element
    Group(T& element)
    {
        Add(element);
    }

    //! @brief Create a group containing multiple elements
    //! @param elements vector of elements
    //! @remark std::reference_wrapper for value semantics
    Group(std::vector<std::reference_wrapper<T>> elements)
    {
        for (auto element : elements)
            Add(element);
    }

    //! @brief Add element to group
    //! @param element single element to add
    void Add(T& element)
    {
        if (mUniqueData.Add(element))
            mData.push_back(&element);
    }

    //! @brief True if element is contained in group
    //! @param element single element to check
    bool Contains(const T& element) const
    {
        return mUniqueData.Contains(element);
    }

    //! @brief True if group is empty
    bool Empty() const
    {
        return Size() == 0;
    }

    //! @brief Number of constituents
    auto Size() const
    {
        return mData.size();
    }

    //! @brief Iterate over group elements (by value not pointer)
    ConstGroupIterator begin() const
    {
        return mData.begin();
    }

    //! @brief Iterate over group elements (by value not pointer)
    ConstGroupIterator end() const
    {
        return mData.end();
    }

    //! @brief Iterate over group elements (by value not pointer)
    GroupIterator begin()
    {
        return mData.begin();
    }

    //! @brief Iterate over group elements (by value not pointer)
    GroupIterator end()
    {
        return mData.end();
    }

private:
    SortedUniqueVector mUniqueData;
    Data mData;
};

//! @brief Unite two groups
template <typename T>
Group<T> Unite(const Group<T>& one, const Group<T>& two)
{
    Group<T> newGroup = one;
    for (auto& t : two)
        newGroup.Add(t);
    return newGroup;
}

//! @brief Returns group with elements of group one that are not in group two
template <typename T>
Group<T> Difference(const Group<T>& one, const Group<T>& two)
{
    Group<T> newGroup;
    for (auto& t : one)
        if (not two.Contains(t))
            newGroup.Add(t);
    return newGroup;
}

//! @brief Returns group with elements that are in both groups
template <typename T>
Group<T> Intersection(const Group<T>& one, const Group<T>& two)
{
    Group<T> newGroup;
    for (auto& t : one)
        if (two.Contains(t))
            newGroup.Add(t);
    return newGroup;
}

//! @brief Returns group with elements that are only in one group not in both
template <typename T>
Group<T> SymmetricDifference(const Group<T>& one, const Group<T>& two)
{
    return Difference(Unite(one, two), Intersection(one, two));
}

} // namespace Group
} // namespace NuTo
