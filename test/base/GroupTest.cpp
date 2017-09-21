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
    group.AddMember(a);
    BOOST_CHECK(group.Contains(a));
    group.AddMember(a);
    BOOST_CHECK(group.Contains(a));
}


struct GroupTestFixture
{
    GroupTestFixture()
    {
        groupOne.AddMember(a);
        groupOne.AddMember(b);
        groupOne.AddMember(c);
        groupOne.AddMember(a); // duplicate!

        groupTwo.AddMember(c);
        groupTwo.AddMember(d);
        groupTwo.AddMember(d); // duplicate!
    }

    Foo a, b, c, d;

    Groups::Group<Foo> groupOne;
    Groups::Group<Foo> groupTwo;
};

BOOST_FIXTURE_TEST_CASE(unionTest, GroupTestFixture)
{
    auto unionGroup = Groups::Unite(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(unionGroup.size(), 4);
    BOOST_CHECK(unionGroup.Contains(a));
    BOOST_CHECK(unionGroup.Contains(b));
    BOOST_CHECK(unionGroup.Contains(c));
    BOOST_CHECK(unionGroup.Contains(d));
}


BOOST_FIXTURE_TEST_CASE(differenceTest, GroupTestFixture)
{
    auto diffGroup = Groups::Difference(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(diffGroup.size(), 2);
    BOOST_CHECK(diffGroup.Contains(a));
    BOOST_CHECK(diffGroup.Contains(b));
    BOOST_CHECK(not diffGroup.Contains(c));
    BOOST_CHECK(not diffGroup.Contains(d));
}


BOOST_FIXTURE_TEST_CASE(intersectionTest, GroupTestFixture)
{
    auto intersectGroup = Groups::Intersection(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(intersectGroup.size(), 1);
    BOOST_CHECK(not intersectGroup.Contains(a));
    BOOST_CHECK(not intersectGroup.Contains(b));
    BOOST_CHECK(intersectGroup.Contains(c));
    BOOST_CHECK(not intersectGroup.Contains(d));
}


BOOST_FIXTURE_TEST_CASE(symmetricDiffTest, GroupTestFixture)
{
    auto symmetricDiffGroup = Groups::SymmetricDifference(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(symmetricDiffGroup.size(), 3);
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

    Groups::Group<int, std::less<int*>> g0;
    g0.AddMember(a);
    g0.AddMember(b);
    g0.AddMember(c);
    Groups::Group<int, std::less<int*>> g1;
    g1.AddMember(b);

    auto intersection = Groups::Intersection(g0, g1);
    BOOST_CHECK(intersection.Contains(b));
}
