#pragma once

#include <algorithm>
#include <vector>
#include <boost/iterator/indirect_iterator.hpp>

namespace NuTo
{
namespace Groups
{

template <typename T>
struct IdCompare
{
    bool operator()(const T* a, const T* b) const
    {
        return a->Id() < b->Id();
    }
};

template <typename T, typename TCompare = IdCompare<T>>
class Group : private std::vector<T*>
{
public:
    typedef std::vector<T*> parent;
    typedef boost::indirect_iterator<typename parent::iterator> GroupIterator;

    void Add(T& element)
    {
        auto it = std::lower_bound(pbegin(), pend(), &element, TCompare());
        if (it == pend() || TCompare()(&element, *it))
            parent::insert(it, &element);
    }

    GroupIterator begin()
    {
        return parent::begin();
    }

    GroupIterator end()
    {
        return parent::end();
    }

    typename parent::iterator pbegin()
    {
        return parent::begin();
    }

    typename parent::iterator pend()
    {
        return parent::end();
    }

    typename parent::const_iterator pcbegin() const
    {
        return parent::cbegin();
    }

    typename parent::const_iterator pcend() const
    {
        return parent::cend();
    }

    bool Contains(const T& element) const
    {
        auto result = std::lower_bound(pcbegin(), pcend(), &element, TCompare());
        if (result == pcend())
            return false;
        return *result == &element;
    }

    using parent::size;
    using typename parent::iterator;
    using typename parent::value_type;
    using parent::insert;
};


template <typename T, typename TCompare>
Group<T, TCompare> Unite(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_union(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator, TCompare());
    return newGroup;
}


template <typename T, typename TCompare>
Group<T, TCompare> Difference(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_difference(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator, TCompare());
    return newGroup;
}


template <typename T, typename TCompare>
Group<T, TCompare> Intersection(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_intersection(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator, TCompare());
    return newGroup;
}


template <typename T, typename TCompare>
Group<T, TCompare> SymmetricDifference(const Group<T, TCompare>& one, const Group<T, TCompare>& two)
{
    Group<T, TCompare> newGroup;
    std::insert_iterator<Group<T, TCompare>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_symmetric_difference(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator,
                                  TCompare());
    return newGroup;
}

} // namespace Group
} // namespace NuTo
