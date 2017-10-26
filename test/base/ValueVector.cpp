#include "BoostUnitTest.h"
#include "base/ValueVector.h"


BOOST_AUTO_TEST_CASE(ValueVectorValidReferences)
{
    NuTo::ValueVector<int> v;
    int& firstElement = v.Add(42);
    for (int i = 0; i < 1e6; ++i)
        v.Add(i);

    int& lastElement = v.Add(4);

    v.Erase(v.begin() + 100, v.end() - 2);

    BOOST_CHECK_EQUAL(v.Size(), 102);

    BOOST_CHECK_EQUAL(firstElement, 42);
    BOOST_CHECK_EQUAL(lastElement, 4);

    v.Erase(std::find(v.begin(), v.end(), firstElement));
    BOOST_CHECK_EQUAL(v.Size(), 101);
}

BOOST_AUTO_TEST_CASE(ValueVectorForwarding)
{
    struct Bar
    {
        Bar(int)
        {
        }
    };

    struct Foo
    {
        Foo(int a, Bar)
            : m(a)
        {
        }
        int m;
    };

    NuTo::ValueVector<Foo> v;
    v.Add(4, Bar(42));
    v.Add(4, 42);

    BOOST_CHECK_EQUAL(v[0].m, 4);
}

BOOST_AUTO_TEST_CASE(ValueVectorRangeLoop)
{
    NuTo::ValueVector<int> v;
    for (int i = 0; i < 10; ++i)
        v.Add(i);

    int sumValue = 0;
    for (int a : v)
        sumValue += a;

    int sumReference = 0;
    for (int& a : v)
        sumReference += a;

    int sumConstReference = 0;
    for (const int& a : v)
        sumConstReference += a;

    BOOST_CHECK_EQUAL(sumValue, 45);
    BOOST_CHECK_EQUAL(sumReference, 45);
    BOOST_CHECK_EQUAL(sumConstReference, 45);
}
