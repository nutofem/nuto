#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/MechanicsException.h"
#include "mechanics/constraints/ConstraintLinearNodeDisplacements2D.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

BOOST_AUTO_TEST_CASE(ConstraintNode2D)
{
    fakeit::Mock<NuTo::NodeBase> node;
    Method(node, GetNum) = 2;
    ConstOverloadedMethod(node, GetDof, int(NuTo::Node::eDof, int)).Using(NuTo::Node::eDof::DISPLACEMENTS, 0) = 1;
    ConstOverloadedMethod(node, GetDof, int(NuTo::Node::eDof, int)).Using(NuTo::Node::eDof::DISPLACEMENTS, 1) = 2;


    NuTo::ConstraintLinearNodeDisplacements2D c(&node.get(), Eigen::Vector2d(0.5, 0.5), 42);

    //c.k
}
