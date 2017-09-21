#include "BoostUnitTest.h"
#include "base/Group.h"
#include "base/UniqueId.h"

#include <iostream>

using namespace NuTo;

struct FooBase : NuTo::UniqueId<FooBase>
{
};
struct Foo : FooBase
{
};

BOOST_AUTO_TEST_CASE(ContainsTest)
{
    Groups::Group<Foo> group;
    Foo a;
    BOOST_CHECK(not group.Contains(a));
    group.Add(a);
    BOOST_CHECK(group.Contains(a));
    group.Add(a);
    BOOST_CHECK(group.Contains(a));
}


struct GroupTestFixture
{
    GroupTestFixture()
    {
        groupOne.Add(a);
        groupOne.Add(b);
        groupOne.Add(c);
        groupOne.Add(a); // duplicate!

        groupTwo.Add(c);
        groupTwo.Add(d);
        groupTwo.Add(d); // duplicate!
    }

    Foo a, b, c, d;

    Groups::Group<Foo> groupOne;
    Groups::Group<Foo> groupTwo;
};

BOOST_FIXTURE_TEST_CASE(unionTest, GroupTestFixture)
{
    auto unionGroup = Groups::Utils::Unite(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(unionGroup.Size(), 4);
    BOOST_CHECK(unionGroup.Contains(a));
    BOOST_CHECK(unionGroup.Contains(b));
    BOOST_CHECK(unionGroup.Contains(c));
    BOOST_CHECK(unionGroup.Contains(d));
}


BOOST_FIXTURE_TEST_CASE(differenceTest, GroupTestFixture)
{
    auto diffGroup = Groups::Utils::Difference(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(diffGroup.Size(), 2);
    BOOST_CHECK(diffGroup.Contains(a));
    BOOST_CHECK(diffGroup.Contains(b));
    BOOST_CHECK(not diffGroup.Contains(c));
    BOOST_CHECK(not diffGroup.Contains(d));
}


BOOST_FIXTURE_TEST_CASE(intersectionTest, GroupTestFixture)
{
    auto intersectGroup = Groups::Utils::Intersection(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(intersectGroup.Size(), 1);
    BOOST_CHECK(not intersectGroup.Contains(a));
    BOOST_CHECK(not intersectGroup.Contains(b));
    BOOST_CHECK(intersectGroup.Contains(c));
    BOOST_CHECK(not intersectGroup.Contains(d));
}


BOOST_FIXTURE_TEST_CASE(symmetricDiffTest, GroupTestFixture)
{
    auto symmetricDiffGroup = Groups::Utils::SymmetricDifference(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(symmetricDiffGroup.Size(), 3);
    BOOST_CHECK(symmetricDiffGroup.Contains(a));
    BOOST_CHECK(symmetricDiffGroup.Contains(b));
    BOOST_CHECK(not symmetricDiffGroup.Contains(c));
    BOOST_CHECK(symmetricDiffGroup.Contains(d));
}

BOOST_AUTO_TEST_CASE(CustomCompare)
{
    int a = 0;
    int b = 6174;
    int c = 42;

    Groups::Group<int, std::less<const int*>> g0;
    g0.Add(a);
    g0.Add(b);
    g0.Add(c);
    Groups::Group<int, std::less<const int*>> g1;
    g1.Add(b);

    auto intersection = Groups::Utils::Intersection(g0, g1);
    BOOST_CHECK(intersection.Contains(b));
}

BOOST_AUTO_TEST_CASE(GroupCtors)
{
    Groups::Group<Foo> empty;
    BOOST_CHECK(empty.Empty());

    Foo a, b, c;
    Groups::Group<Foo> one(a);
    BOOST_CHECK(one.Contains(a));

    Groups::Group<Foo> many({a, b, c, a});
    BOOST_CHECK(many.Contains(a));
    BOOST_CHECK(many.Contains(b));
    BOOST_CHECK(many.Contains(c));
    BOOST_CHECK_EQUAL(many.Size(), 3);
}
