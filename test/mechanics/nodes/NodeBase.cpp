#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/MechanicsEnums.h"

using namespace NuTo;
using namespace fakeit;

BOOST_AUTO_TEST_CASE(NodeBaseSetGet)
{
    fakeit::Mock<NodeBase> mockNode;
    NuTo::NodeBase& nd = mockNode.get();
    When(OverloadedMethod(mockNode, Set, void(Node::eDof, int, const Eigen::VectorXd&))).AlwaysReturn();

    nd.Set(NuTo::Node::eDof::DISPLACEMENTS, 2.7);
    Verify(OverloadedMethod(mockNode, Set, void(Node::eDof, int, const Eigen::VectorXd&))
                   .Matching([](auto, auto timeDerivative, const auto&) { return timeDerivative == 0; }))
            .Once();

    nd.Set(NuTo::Node::eDof::DISPLACEMENTS, 1, 5.6);
    Verify(OverloadedMethod(mockNode, Set, void(Node::eDof, int, const Eigen::VectorXd&))
                   .Matching([](auto, auto timeDerivative, const auto&) { return timeDerivative == 1; }))
            .Once();
}
