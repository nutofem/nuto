#include "BoostUnitTest.h"
#include "math/Average.h"
#include <random>

constexpr long double e = 1.e-12;

BOOST_AUTO_TEST_CASE(AverageBasics)
{
    NuTo::Average a;
    BOOST_CHECK_EQUAL(a.mNum, 0);
    BOOST_CHECK_SMALL(a.mCurrentAverage, e);

    a(5.);
    BOOST_CHECK_EQUAL(a.mNum, 1);
    BOOST_CHECK_CLOSE(a.mCurrentAverage, 5., e);

    a(5.);
    BOOST_CHECK_EQUAL(a.mNum, 2);
    BOOST_CHECK_CLOSE(a.mCurrentAverage, 5., e);

    a(2.);
    BOOST_CHECK_EQUAL(a.mNum, 3);
    BOOST_CHECK_CLOSE(a.mCurrentAverage, 4., e);
}

BOOST_AUTO_TEST_CASE(AveragePrecision)
{
    NuTo::Average a;
    constexpr auto n = 1e6;
    std::mt19937 rng;
    std::uniform_int_distribution<int> ints(0, n);
    long int sum = 0;
    for (auto i = 0; i < n; ++i)
    {
        auto number = ints(rng);
        sum += number;
        a(number);
    }
    long double mean = static_cast<double>(sum) / n;
    BOOST_CHECK_CLOSE_FRACTION(a.mCurrentAverage, mean, e);
}
