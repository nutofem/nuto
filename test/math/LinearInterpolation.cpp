#include "BoostUnitTest.h"

#include <array>
#include "nuto/math/LinearInterpolation.h"
#include "nuto/base/Exception.h"


#define CHECK_CLOSE_COLLECTION(aa, bb, tolerance)                                                                      \
    {                                                                                                                  \
        using std::distance;                                                                                           \
        using std::begin;                                                                                              \
        using std::end;                                                                                                \
        auto a = begin(aa), ae = end(aa);                                                                              \
        auto b = begin(bb);                                                                                            \
        BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(bb)));                                                    \
        for (; a != ae; ++a, ++b)                                                                                      \
        {                                                                                                              \
            BOOST_CHECK_CLOSE(*a, *b, tolerance);                                                                      \
        }                                                                                                              \
    }

BOOST_AUTO_TEST_CASE(linear_interpolation)
{
    std::vector<std::array<double, 2>> values = {{{0.0, 0.0}, {100.0, 1.0}}};
    auto interpolation = NuTo::Math::LinearInterpolation(values);

    BOOST_CHECK_CLOSE(interpolation(50.0), 0.5, 1e-10);
    BOOST_CHECK_CLOSE(interpolation(0.0), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(interpolation(100.0), 1.0, 1e-10);

    double derivative = interpolation.derivative(50.0);
    BOOST_CHECK_CLOSE(derivative, 1e-2, 1e-10);

    BOOST_CHECK_THROW(interpolation(-1.0), NuTo::Exception);
    BOOST_CHECK_THROW(interpolation(101.0), NuTo::Exception);

    values = {{{0.0, 0.0}}};
    BOOST_CHECK_THROW(auto interpolation = NuTo::Math::LinearInterpolation(values), NuTo::Exception);

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
