#pragma once

#include <boost/range/any_range.hpp>

namespace NuTo
{
template <typename T>
using Base = boost::any_range<T, boost::forward_traversal_tag, T&, std::ptrdiff_t>;

//! Provides a view to a collection independent from its underlying container using polymorphic iterators provided by
//! boost::any_range
template <typename T>
class ContainerView : public Base<T>
{
public:
    ContainerView(const ContainerView&) = default;


    ContainerView(Base<T> other)
        : Base<T>(other)
    {
    }

    template <typename TOther>
    ContainerView(const TOther& other)
        : Base<T>(other.begin(), other.end())
    {
    }

    ContainerView(std::initializer_list<std::reference_wrapper<T>> list)
        : Base<T>(list.begin(), list.end())
        , m(std::move(list))
    {
    }

private:
    std::initializer_list<std::reference_wrapper<T>> m;
};
} /* NuTo */
