#pragma once

#include <boost/ptr_container/ptr_vector.hpp>

namespace NuTo
{

//! @brief container that stores values of T and keeps references to these values valid
template <typename T>
class ValueVector : public boost::ptr_vector<T>
{
public:

    template <typename... TArgs>
    T& Add(TArgs&&... args)
    {
        this->push_back(new T(std::forward<TArgs>(args)...));
        return this->back();
    }

    T& Add(T&& t)
    {
        this->push_back(new T(t));
        return this->back();
    }
};
} /* NuTo */
