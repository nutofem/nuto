#include "BoostUnitTest.h"
#include "base/ValueVector.h"


BOOST_AUTO_TEST_CASE(ValueVectorValidReferences)
{
    NuTo::ValueVector<int> v;
    int& firstElement = v.Add(42);
    for (int i = 0; i < 1e6; ++i)
        v.Add(i);

    int& lastElement = v.Add(4);
    v.erase_if([](int a) { return a >= 100; });
    BOOST_CHECK_EQUAL(v.size(), 102); // 100 remaining elements + first(42) + last(4)

    BOOST_CHECK_EQUAL(firstElement, 42);
    BOOST_CHECK_EQUAL(lastElement, 4);
}

BOOST_AUTO_TEST_CASE(ValueVectorForwarding)
{
    struct Bar
    {
        Bar(int) {}
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

    BOOST_CHECK_EQUAL(v.front().m, 4); 
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
