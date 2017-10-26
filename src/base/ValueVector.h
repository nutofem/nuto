#pragma once

#include <vector>
#include <memory>
#include <boost/iterator/indirect_iterator.hpp>


namespace NuTo
{

//! @brief container that stores values of T and keeps references to these values valid
//! @remark Normally, it is a no-brainer to use boost::ptr_vector here. It provides all
//! the nice value semantics for free. However, there is an issue regarding move.
//!         https://stackoverflow.com/questions/44130076/moving-a-boostptr-vector
//! Therefore, here is a custom implementation with only the basic features.
template <typename T>
class ValueVector
{
    using Data = std::vector<std::unique_ptr<T>>;
    //! @brief indirect (dereferencing) iterator to provide value semantics for the iterators
    typedef boost::indirect_iterator<typename Data::iterator> IndirectIterator;
    typedef boost::indirect_iterator<typename Data::const_iterator> ConstIndirectIterator;

public:
    template <typename... TArgs>
    T& Add(TArgs&&... args)
    {
        mData.push_back(std::make_unique<T>(std::forward<TArgs>(args)...));
        return *mData.back();
    }

    T& Add(T&& t)
    {
        mData.push_back(std::make_unique<T>(std::move(t)));
        return *mData.back();
    }

    IndirectIterator begin()
    {
        return mData.begin();
    }

    IndirectIterator end()
    {
        return mData.end();
    }

    ConstIndirectIterator begin() const
    {
        return mData.begin();
    }

    ConstIndirectIterator end() const
    {
        return mData.end();
    }

    typename Data::size_type Size() const
    {
        return mData.size();
    }

    const T& operator[](int i) const
    {
        return *mData[i];
    }
    T& operator[](int i)
    {
        return *mData[i];
    }

    IndirectIterator Erase(IndirectIterator it)
    {
        return mData.erase(it.base());
    }

    IndirectIterator Erase(IndirectIterator from, IndirectIterator to)
    {
        return mData.erase(from.base(), to.base());
    }

private:
    Data mData;
};
} /* NuTo */
