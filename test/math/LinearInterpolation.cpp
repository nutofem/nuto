#define BOOST_TEST_MODULE LinearInterpolation
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <array>
#include "math/LinearInterpolation.h"
#include "math/MathException.h"

// needed for building with clang when boost test has been built with gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size()));
}

#define CHECK_CLOSE_COLLECTION(aa, bb, tolerance) { \
    using std::distance; \
    using std::begin; \
    using std::end; \
    auto a = begin(aa), ae = end(aa); \
    auto b = begin(bb); \
    BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(bb))); \
    for(; a != ae; ++a, ++b) { \
        BOOST_CHECK_CLOSE(*a, *b, tolerance); \
    } \
}

BOOST_AUTO_TEST_CASE(linear_interpolation)
{
    std::vector<std::array<double, 2>> values = {{{0.0, 0.0},
                                                {100.0, 1.0}}};
    auto interpolation = NuTo::Math::LinearInterpolation(values);

    BOOST_CHECK_CLOSE(interpolation(50.0), 0.5, 1e-10);
    BOOST_CHECK_CLOSE(interpolation(0.0), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(interpolation(100.0), 1.0, 1e-10);

    double derivative = interpolation.derivative(50.0);
    BOOST_CHECK_CLOSE(derivative, 1e-2, 1e-10);

    BOOST_CHECK_THROW(interpolation(-1.0), NuTo::out_of_range);
    BOOST_CHECK_THROW(interpolation(101.0), NuTo::out_of_range);

    values = {{{0.0, 0.0}}};
    BOOST_CHECK_THROW(auto interpolation = NuTo::Math::LinearInterpolation(values), NuTo::invalid_argument);

    // values need not be in the correct order
    values = {{{100.0, 1.0}, {0.0, 0.0}}};
    BOOST_CHECK_CLOSE(interpolation(50.0), 0.5, 1e-10);
    BOOST_CHECK_CLOSE(interpolation(0.0), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(interpolation(100.0), 1.0, 1e-10);

    // function object
    std::function<double(double)> f = interpolation.f;
    BOOST_CHECK_CLOSE(f(50.0), 0.5, 1e-10);
    BOOST_CHECK_CLOSE(f(0.0), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(f(100.0), 1.0, 1e-10);

    std::function<double(double)> f_prime = interpolation.df;
    BOOST_CHECK_CLOSE(f_prime(50.0), 1e-2, 1e-10);
}
