#include "BoostUnitTest.h"
#include "base/Group.h"

#include <iostream>

using namespace NuTo;


BOOST_AUTO_TEST_CASE(ContainsTest)
{
    Groups::Group<int> group;
    int a = 0;
    BOOST_CHECK(not group.Contains(a));
    group.AddMember(a);
    BOOST_CHECK(group.Contains(a));
    group.AddMember(a);
    BOOST_CHECK(group.Contains(a));
}

class Foo
{
};

struct GroupTestFixture
{
    GroupTestFixture()
    {
        groupOne.AddMember(a);
        groupOne.AddMember(b);
        groupOne.AddMember(c);

        groupTwo.AddMember(c);
        groupTwo.AddMember(d);
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
}


BOOST_FIXTURE_TEST_CASE(intersectionTest, GroupTestFixture)
{
    auto intersectGroup = Groups::Intersection(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(intersectGroup.size(), 1);
    BOOST_CHECK(intersectGroup.Contains(c));
}


BOOST_FIXTURE_TEST_CASE(symmetricDiffTest, GroupTestFixture)
{
    auto symmetricDiffGroup = Groups::SymmetricDifference(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(symmetricDiffGroup.size(), 3);
    BOOST_CHECK(symmetricDiffGroup.Contains(a));
    BOOST_CHECK(symmetricDiffGroup.Contains(b));
    BOOST_CHECK(symmetricDiffGroup.Contains(d));
}
