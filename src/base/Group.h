#pragma once

#include <algorithm>
#include <set>
#include <boost/iterator/indirect_iterator.hpp>

namespace NuTo
{
namespace Groups
{

template <typename T>
class Group : private std::set<T*>
{
public:
    typedef std::set<T*> parent;
    typedef boost::indirect_iterator<typename parent::iterator> GroupIterator;

    void AddMember(T& element)
    {
        parent::insert(&element);
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

    bool Contains(const T& element)
    {
        auto result = std::find(parent::begin(), parent::end(), &element);
        if (result != parent::end())
            return true;
        else
            return false;
    }

    using parent::size;
    using typename parent::iterator;
    using typename parent::value_type;
    using parent::insert;
};


template <typename T>
Group<T> Unite(const Group<T>& one, const Group<T>& two)
{
    Group<T> newGroup;
    std::insert_iterator<Group<T>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_union(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator);
    return newGroup;
}


template <typename T>
Group<T> Difference(const Group<T>& one, const Group<T>& two)
{
    Group<T> newGroup;
    std::insert_iterator<Group<T>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_difference(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator);
    return newGroup;
}


template <typename T>
Group<T> Intersection(const Group<T>& one, const Group<T>& two)
{
    Group<T> newGroup;
    std::insert_iterator<Group<T>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_intersection(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator);
    return newGroup;
}


template <typename T>
Group<T> SymmetricDifference(const Group<T>& one, const Group<T>& two)
{
    Group<T> newGroup;
    std::insert_iterator<Group<T>> newGroupInsertIterator(newGroup, newGroup.pbegin());
    std::set_symmetric_difference(one.pcbegin(), one.pcend(), two.pcbegin(), two.pcend(), newGroupInsertIterator);
    return newGroup;
}

} // namespace Group
} // namespace NuTo
