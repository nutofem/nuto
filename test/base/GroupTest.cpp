#include "BoostUnitTest.h"
#include "base/Group.h"

using namespace NuTo::Groups;

struct Foo
{
};

bool operator==(const Foo& a, const Foo& b)
{
    return &a == &b;
}

BOOST_AUTO_TEST_CASE(ContainsTest)
{
    Group<Foo> group;
    Foo a, b, c, d;
    BOOST_CHECK(not group.Contains(a));
    group.Add(a);
    BOOST_CHECK(group.Contains(a));
    group.Add(a);
    BOOST_CHECK(group.Contains(a));
    group.Add(b);
    group.Add(c);
    group.Add(d);
    BOOST_CHECK(group.Contains(c));
}


struct GroupTestFixture
{
    GroupTestFixture()
    {
        groupOne.Add(a);
        groupOne.Add(b);
        groupOne.Add(c);
        groupOne.Add(b);
        groupOne.Add(a);

        groupTwo.Add(e);
        groupTwo.Add(d);
        groupTwo.Add(d);
        groupTwo.Add(c);
        groupTwo.Add(d);
        groupTwo.Add(e);
        groupTwo.Add(e);
        groupTwo.Add(c);
    }

    Foo a, b, c, d, e;

    Group<Foo> groupOne;
    Group<Foo> groupTwo;
};

BOOST_FIXTURE_TEST_CASE(unionTest, GroupTestFixture)
{
    auto unionGroup = Unite(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(unionGroup.Size(), 5);
    auto it = unionGroup.begin();
    BOOST_CHECK(*it == a);

    std::advance(it, 1);
    BOOST_CHECK(*it == b);

    std::advance(it, 1);
    BOOST_CHECK(*it == c);

    std::advance(it, 1);
    BOOST_CHECK(*it == e);

    std::advance(it, 1);
    BOOST_CHECK(*it == d);
}


BOOST_FIXTURE_TEST_CASE(differenceTest, GroupTestFixture)
{
    auto diffGroup = Difference(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(diffGroup.Size(), 2);
    BOOST_CHECK(diffGroup.Contains(a));
    BOOST_CHECK(diffGroup.Contains(b));
    BOOST_CHECK(not diffGroup.Contains(c));
    BOOST_CHECK(not diffGroup.Contains(d));

    auto it = diffGroup.begin();
    BOOST_CHECK(*it == a);
    BOOST_CHECK(*std::next(it) == b);
}


BOOST_FIXTURE_TEST_CASE(intersectionTest, GroupTestFixture)
{
    auto intersectGroup = Intersection(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(intersectGroup.Size(), 1);
    BOOST_CHECK(not intersectGroup.Contains(a));
    BOOST_CHECK(not intersectGroup.Contains(b));
    BOOST_CHECK(intersectGroup.Contains(c));
    BOOST_CHECK(not intersectGroup.Contains(d));
}


BOOST_FIXTURE_TEST_CASE(symmetricDiffTest, GroupTestFixture)
{
    auto symmetricDiffGroup = SymmetricDifference(groupOne, groupTwo);
    BOOST_CHECK_EQUAL(symmetricDiffGroup.Size(), 4);
    BOOST_CHECK(symmetricDiffGroup.Contains(a));
    BOOST_CHECK(symmetricDiffGroup.Contains(b));
    BOOST_CHECK(not symmetricDiffGroup.Contains(c));
    BOOST_CHECK(symmetricDiffGroup.Contains(d));
    BOOST_CHECK(symmetricDiffGroup.Contains(e));
}

BOOST_AUTO_TEST_CASE(CustomCompare)
{
    int a = 0;
    int b = 6174;
    int c = 42;

    Group<int> g0;
    g0.Add(a);
    g0.Add(b);
    g0.Add(c);
    Group<int> g1;
    g1.Add(b);

    auto intersection = Intersection(g0, g1);
    BOOST_CHECK(intersection.Contains(b));
}

BOOST_AUTO_TEST_CASE(GroupCtors)
{
    Group<Foo> empty;
    BOOST_CHECK(empty.Empty());

    Foo a, b, c;
    Group<Foo> one(a);
    BOOST_CHECK(one.Contains(a));

    Group<Foo> many({a, b, c, a});
    BOOST_CHECK(many.Contains(a));
    BOOST_CHECK(many.Contains(b));
    BOOST_CHECK(many.Contains(c));
    BOOST_CHECK_EQUAL(many.Size(), 3);
}

void BraceInitialization(const Group<Foo>& group)
{
    BOOST_CHECK_EQUAL(group.Size(), 2); 
}

BOOST_AUTO_TEST_CASE(GroupBraceInitialize)
{
    Foo a, b;
    BraceInitialization({a, b});
}
