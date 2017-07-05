#include "BoostUnitTest.h"

#include <array>
#include "math/CubicSplineInterpolation.h"
#include "math/MathException.h"

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

BOOST_AUTO_TEST_CASE(cubic_spline_interpolation)
{
    std::vector<std::array<double, 2>> values = {{{0.0, 0.0}, {100.0, 1.0}}};
    auto interpolation = NuTo::Math::CubicSplineInterpolation(values);

    BOOST_CHECK_CLOSE(interpolation(50.0), 0.5, 1e-10);

    std::vector<double> x = {1.0, 3.0, 5.0, 7.0, 9.0};
    std::vector<double> y = {1.0, 3.0, 5.0, 7.0, 9.0};
    std::vector<std::array<double, 2>> values2;
    values2.resize(x.size());
    for (unsigned i = 0; i < x.size(); ++i)
    {
        values2[i][0] = x[i];
        values2[i][1] = y[i];
    }

    auto interpolation2 = NuTo::Math::CubicSplineInterpolation(values2);
    std::vector<double> output;
    output.resize(x.size());
    for (unsigned i = 0; i < x.size(); ++i)
    {
        output[i] = interpolation2(x[i]);
    }
    CHECK_CLOSE_COLLECTION(output, y, 1e-10);
}
