#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"

#define BOOST_TEST_MODULE ContinuumElementTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <fakeit.hpp>

using namespace NuTo;
using namespace fakeit;

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

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
