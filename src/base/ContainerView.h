#pragma once

#include <boost/range/any_range.hpp>

namespace NuTo
{
//! Provides a view to a collection independent from its underlying container using polymorphic iterators provided by
//! boost::any_range
template <typename T>
using ContainerView = boost::any_range<T, boost::forward_traversal_tag, T&, std::ptrdiff_t>;
} /* NuTo */
