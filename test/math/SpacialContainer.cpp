//
// Created by Thomas Titscher on 2/8/17.
//

#include "BoostUnitTest.h"
#include "math/SpacialContainer.h"


std::vector<int> GetValues()
{
    std::vector<int> vals;
    vals.push_back(0);
    vals.push_back(11);
    vals.push_back(20);
    vals.push_back(1101);
    return vals;
}

std::vector<Eigen::Vector2d> GetCoords()
{
    std::vector<Eigen::Vector2d> coords;
    coords.push_back(Eigen::Vector2d({0,0}));
    coords.push_back(Eigen::Vector2d({1,1}));
    coords.push_back(Eigen::Vector2d({2,0}));
    coords.push_back(Eigen::Vector2d({1,1.01}));
    return coords;
}

struct TestContainer
{
    TestContainer() : mValues(GetValues()), s(GetCoords(), mValues) {}
    std::vector<int> mValues;
    NuTo::SpacialContainer<int, 2> s;
};


BOOST_AUTO_TEST_CASE(SpacialContainer_GetDuplicateIDs)
{
    TestContainer test;
    auto ids0 = test.s.GetDuplicateIDs(0, 0.1);
    BOOST_CHECK_EQUAL(ids0.size(), 1);
    BOOST_CHECK_EQUAL(ids0[0], 0);

    auto ids1 = test.s.GetDuplicateIDs(1, 0.1);
    BOOST_CHECK_EQUAL(ids1.size(), 2);
    BOOST_CHECK_EQUAL(ids1[0], 1);
    BOOST_CHECK_EQUAL(ids1[1], 3);
}

BOOST_AUTO_TEST_CASE(SpacialContainer_GetAllDuplicates)
{
    TestContainer test;
    auto ids = test.s.GetAllDuplicateIDs(0.1);
    
    BOOST_CHECK_EQUAL(ids.size(), 3);

    BOOST_CHECK_EQUAL(ids[0].size(), 1);
    BOOST_CHECK_EQUAL(ids[0][0], 0);

    BOOST_CHECK_EQUAL(ids[1].size(), 2);
    BOOST_CHECK_EQUAL(ids[1][0], 1);
    BOOST_CHECK_EQUAL(ids[1][1], 3);

    BOOST_CHECK_EQUAL(ids[2].size(), 1);
    BOOST_CHECK_EQUAL(ids[2][0], 2);
}

BOOST_AUTO_TEST_CASE(SpacialContainer_GetDuplicatesValues)
{
    TestContainer test;
    auto duplicates = test.s.GetAllDuplicateValues(0.1);

    BOOST_CHECK_EQUAL(duplicates.size(), 3);

    BOOST_CHECK_EQUAL(duplicates[0].size(), 1);
    BOOST_CHECK_EQUAL(duplicates[0][0], 0);

    BOOST_CHECK_EQUAL(duplicates[1].size(), 2);
    BOOST_CHECK_EQUAL(duplicates[1][0], 11);
    BOOST_CHECK_EQUAL(duplicates[1][1], 1101);

    BOOST_CHECK_EQUAL(duplicates[2].size(), 1);
    BOOST_CHECK_EQUAL(duplicates[2][0], 20);
}

