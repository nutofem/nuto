//
// Created by Thomas Titscher on 2/8/17.
//

#include "BoostUnitTest.h"
#include "math/SpatialContainer.h"

struct TestNode
{
    int value;
    Eigen::VectorXd coord;
};

struct TestNodeCoord
{
    Eigen::VectorXd operator()(const TestNode& r)
    {
        return r.coord;
    }
};

std::vector<TestNode> GetTestNodes()
{
    std::vector<TestNode> v;
    v.push_back({0, Eigen::Vector2d(0, 0)});
    v.push_back({11, Eigen::Vector2d(1, 1)});
    v.push_back({20, Eigen::Vector2d(2, 0)});
    v.push_back({1101, Eigen::Vector2d(1, 1.01)});
    return v;
}


BOOST_AUTO_TEST_CASE(SpatialContainer_GetDuplicateIDs)
{
    auto testNodes = GetTestNodes();
    NuTo::SpatialContainer<TestNode, TestNodeCoord> s(testNodes);

    auto ids0 = s.FindIDsWithinRadius(0, 0.1);
    BOOST_CHECK_EQUAL(ids0.size(), 1);
    BOOST_CHECK_EQUAL(ids0[0], 0);

    auto ids1 = s.FindIDsWithinRadius(1, 0.1);
    BOOST_CHECK_EQUAL(ids1.size(), 2);
    BOOST_CHECK_EQUAL(ids1[0], 1);
    BOOST_CHECK_EQUAL(ids1[1], 3);
}

BOOST_AUTO_TEST_CASE(SpatialContainer_GetAllDuplicates)
{
    auto testNodes = GetTestNodes();
    NuTo::SpatialContainer<TestNode, TestNodeCoord> s(testNodes);

    auto ids = s.GetAllDuplicateIDs(0.1);

    BOOST_CHECK_EQUAL(ids.size(), 3);

    BOOST_CHECK_EQUAL(ids[0].size(), 1);
    BOOST_CHECK_EQUAL(ids[0][0], 0);

    BOOST_CHECK_EQUAL(ids[1].size(), 2);
    BOOST_CHECK_EQUAL(ids[1][0], 1);
    BOOST_CHECK_EQUAL(ids[1][1], 3);

    BOOST_CHECK_EQUAL(ids[2].size(), 1);
    BOOST_CHECK_EQUAL(ids[2][0], 2);
}

BOOST_AUTO_TEST_CASE(SpatialContainer_GetDuplicatesValues)
{
    auto testNodes = GetTestNodes();
    NuTo::SpatialContainer<TestNode, TestNodeCoord> s(testNodes);

    auto duplicates = s.GetAllDuplicateValues(0.1);

    BOOST_CHECK_EQUAL(duplicates.size(), 3);

    BOOST_CHECK_EQUAL(duplicates[0].size(), 1);
    BOOST_CHECK_EQUAL(duplicates[0][0].value, 0);

    BOOST_CHECK_EQUAL(duplicates[1].size(), 2);
    BOOST_CHECK_EQUAL(duplicates[1][0].value, 11);
    BOOST_CHECK_EQUAL(duplicates[1][1].value, 1101);

    BOOST_CHECK_EQUAL(duplicates[2].size(), 1);
    BOOST_CHECK_EQUAL(duplicates[2][0].value, 20);
}

BOOST_AUTO_TEST_CASE(SpatialContainer_HasEntryAtCoordinate)
{
    auto testNodes = GetTestNodes();
    NuTo::SpatialContainer<TestNode, TestNodeCoord> s(testNodes);

    BOOST_CHECK(s.HasEntryAtCoordinate(Eigen::Vector2d({0, 0}), 0.1));
    BOOST_CHECK(not s.HasEntryAtCoordinate(Eigen::Vector2d({1, 0}), 0.1));
}