#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "math/SparseMatrixCSRVector2General.h"

const NuTo::Node::eDof eDofDisp = NuTo::Node::eDof::DISPLACEMENTS;

auto GetMockNode()
{
    fakeit::Mock<NuTo::NodeBase> node;
    Method(node, GetNum) = 2;
    Method(node, GetNumDofs) = 2;
    ConstOverloadedMethod(node, GetDof, int(NuTo::Node::eDof, int)).Using(eDofDisp, 0) = 1;
    ConstOverloadedMethod(node, GetDof, int(NuTo::Node::eDof, int)).Using(eDofDisp, 1) = 2;
    return node;
}

BOOST_AUTO_TEST_CASE(ConstraintX)
{
    fakeit::Mock<NuTo::NodeBase> node = GetMockNode();

    NuTo::Constraint::Constraints c;
    c.Add(eDofDisp, NuTo::Constraint::Component(node.get(), {NuTo::eDirection::X}));
    BOOST_CHECK_EQUAL(c.GetNumEquations(eDofDisp), 1);
    
    NuTo::SparseMatrixCSRVector2General<double> m(1, 3);
    c.BuildConstraintMatrix(m, eDofDisp);
    BoostUnitTest::CheckEigenMatrix(m.ConvertToFullMatrix(), Eigen::Vector3d(0,1,0).transpose());
}

BOOST_AUTO_TEST_CASE(ConstraintXandY)
{
    fakeit::Mock<NuTo::NodeBase> node = GetMockNode();

    NuTo::Constraint::Constraints c;
    c.Add(eDofDisp, NuTo::Constraint::Component(node.get(), {NuTo::eDirection::X, NuTo::eDirection::Y}));
    BOOST_CHECK_EQUAL(c.GetNumEquations(eDofDisp), 2);

    NuTo::SparseMatrixCSRVector2General<double> m(2, 3);
    c.BuildConstraintMatrix(m, eDofDisp);
    Eigen::MatrixXd expected(2,3);
    expected.setZero();
    expected(0,1) = 1;
    expected(1,2) = 1;
    BoostUnitTest::CheckEigenMatrix(m.ConvertToFullMatrix(), expected); 
    
    std::cout << c << std::endl;
}

BOOST_AUTO_TEST_CASE(ConstraintDirection)
{
    fakeit::Mock<NuTo::NodeBase> node = GetMockNode();

    NuTo::Constraint::Constraints c;
    c.Add(eDofDisp, NuTo::Constraint::Direction(node.get(), Eigen::Vector2d(1,-1), NuTo::Constraint::RhsConstant(0)));
    BOOST_CHECK_EQUAL(c.GetNumEquations(eDofDisp), 1);

    NuTo::SparseMatrixCSRVector2General<double> m(1, 3);
    c.BuildConstraintMatrix(m, eDofDisp);
    Eigen::Vector3d expected(0, 1/sqrt(2), - 1/sqrt(2));
    BoostUnitTest::CheckEigenMatrix(m.ConvertToFullMatrix(), expected.transpose()); 

    std::cout << c << std::endl;
}
