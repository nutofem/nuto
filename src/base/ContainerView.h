#pragma once

#include <boost/range/any_range.hpp>

namespace NuTo
{
template <typename T>
using Base = boost::any_range<T, boost::forward_traversal_tag, T&, std::ptrdiff_t>;

//! Provides a view to a collection independent from its underlying container using polymorphic iterators provided by
//! boost::any_range. This `range` basically wraps `.begin()` and `.end()` of a container.
//! @tparam T data type
template <typename T>
class ContainerView : public Base<T>
{
public:
    //! Copy ctor
    ContainerView(const ContainerView&) = default;

    //! Contstruction from any container that provides `.begin()` and `.end()`
    //! @tparam TOther any container type, e.g. from stl or NuTo::Group
    template <typename TOther>
    ContainerView(TOther&& other)
        : Base<T>(other.begin(), other.end())
    {
    }

    //! Construction from an initializer list
    ContainerView(std::initializer_list<std::reference_wrapper<T>>&& l)
        : Base<T>(l)
    {
    }
};
} /* NuTo */
