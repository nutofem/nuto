#pragma once

#include <algorithm>
#include <vector>
#include <boost/iterator/indirect_iterator.hpp>

namespace NuTo
{
namespace Groups
{

template <typename T>
struct IdSmaller
{
    bool operator()(const T* a, const T* b) const
    {
        return a->Id() < b->Id();
    }
};

class Utils;

//! @brief Ordered container class for elements, nodes and the like. No duplicates.
//!
//! The group saves pointers to its members but can be iterated through
//! by value. By default the members are sorted by id (this must be unique)
//! that is accessed by calling Id(), e.g as provided by NuTo::UniqueId. Other
//! methods operating on the group also use this comparator.
//! This behaviour can be changed by providing a different comparator but the basis
//! of the comparison has to be a unique attribute of the group member.
//! See test file (CutomCompare).
//!
//! @tparam T type of group members
//! @tparam TCompare functor representing the comparison operation used for sorting
template <typename T, typename TCompare = IdSmaller<T>>
class Group : private std::vector<T*>
{
    friend class std::insert_iterator<Group>;

public:
    typedef std::vector<T*> parent;
    typedef boost::indirect_iterator<typename parent::iterator> GroupIterator;

    //! @brief Empty group
    Group() = default;

    //! @brief Group containing one element
    Group(T& element)
    {
        Add(element);
    }

    //! @brief Group containing multiple elements
    Group(std::vector<std::reference_wrapper<T>> elements)
    {
        for (auto element : elements)
            Add(element);
    }

    //! @brief Add element to group
    void Add(T& element)
    {
        auto it = std::lower_bound(pcbegin(), pcend(), &element, TCompare());
        if (it == pcend() || TCompare()(&element, *it))
            parent::insert(it, &element);
    }

    //! @brief True if element is contained in group
    //!
    //! Comment:
    //! Search for element is done using Comparator,
    //! Equality is checked by pointer equivalence.
    bool Contains(const T& element) const
    {
        auto result = std::lower_bound(pcbegin(), pcend(), &element, TCompare());
        if (result == pcend())
            return false;
        return *result == &element;
    }

    //! @brief True if group is empty
    bool Empty() const
    {
        return Size() == 0;
    }

    //! @brief Number of constituents
    auto Size() const
    {
        return size();
    }

    //! @brief Iterate over group elements (by value not pointer)
    GroupIterator begin()
    {
        return parent::begin();
    }

    //! @brief Iterate over group elements (by value not pointer)
    GroupIterator end()
    {
        return parent::end();
    }

    typename parent::iterator pbegin()
    {
        return parent::begin();
    }

    typename parent::const_iterator pcbegin() const
    {
        return parent::cbegin();
    }

    typename parent::const_iterator pcend() const
    {
        return parent::cend();
    }

private:
    using parent::size;
};

//! @brief Unite two groups
template <typename T, typename TCompare>
static Group<T, TCompare> Unite(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_union(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator, TCompare());
    return newGroup;
}

//! @brief Returns group with elements of group one that are not in group two
template <typename T, typename TCompare>
static Group<T, TCompare> Difference(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_difference(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator, TCompare());
    return newGroup;
}

//! @brief Returns group with elements that are in both groups
template <typename T, typename TCompare>
static Group<T, TCompare> Intersection(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_intersection(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator, TCompare());
    return newGroup;
}

//! @brief Returns group with elements that are only in one group not in both
template <typename T, typename TCompare>
static Group<T, TCompare> SymmetricDifference(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_symmetric_difference(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator,
                                  TCompare());
    return newGroup;
}

} // namespace Group
} // namespace NuTo
