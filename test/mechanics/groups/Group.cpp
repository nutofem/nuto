#include "BoostUnitTest.h"
#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"

#include <fakeit.hpp>

using namespace NuTo;
using namespace fakeit;

struct GroupTestFixture
{
    GroupTestFixture()
    {

        Mock<NodeBase> mockNode;

        NodeBase& node1 = mockNode.get();
        NodeBase& node2 = mockNode.get();
        NodeBase& node3 = mockNode.get();
        NodeBase& node4 = mockNode.get();

        nodeGroupOne.AddMember(1, &node1);
        nodeGroupOne.AddMember(2, &node2);
        nodeGroupOne.AddMember(3, &node3);

        nodeGroupTwo.AddMember(3, &node3);
        nodeGroupTwo.AddMember(4, &node4);

    }

    Group<NodeBase> nodeGroupOne;
    Group<NodeBase> nodeGroupTwo;
};


BOOST_FIXTURE_TEST_CASE(unionTest, GroupTestFixture)
{
    auto unionGroup = nodeGroupOne.Unite(&nodeGroupTwo);
    BOOST_CHECK_EQUAL(unionGroup->GetNumMembers(), 4);
}


BOOST_FIXTURE_TEST_CASE(differenceTest, GroupTestFixture)
{
    auto diffGroup = nodeGroupOne.Difference(&nodeGroupTwo);
    BOOST_CHECK_EQUAL(diffGroup->GetNumMembers(), 2);
}


BOOST_FIXTURE_TEST_CASE(intersectionTest, GroupTestFixture)
{
    auto intersectGroup = nodeGroupOne.Intersection(&nodeGroupTwo);
    BOOST_CHECK_EQUAL(intersectGroup->GetNumMembers(), 1);
}


BOOST_FIXTURE_TEST_CASE(symmetricDiffTest, GroupTestFixture)
{
    auto symmetricDiffGroup = nodeGroupOne.SymmetricDifference(&nodeGroupTwo);
    BOOST_CHECK_EQUAL(symmetricDiffGroup->GetNumMembers(), 3);
}
