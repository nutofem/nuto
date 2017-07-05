#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "base/Exception.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constraints/Constraints.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/nodes/NodeBase.h"

const NuTo::Node::eDof eDofDisp = NuTo::Node::eDof::DISPLACEMENTS;

auto GetMockNode(std::vector<int> dofs)
{
    fakeit::Mock<NuTo::NodeBase> node;
    Method(node, GetNum) = dofs.size();
    Method(node, GetNumDofs) = dofs.size();
    ConstOverloadedMethod(node, GetDof, int(NuTo::Node::eDof, int)).Using(eDofDisp, 0) = dofs[0];
    ConstOverloadedMethod(node, GetDof, int(NuTo::Node::eDof, int)).Using(eDofDisp, 1) = dofs[1];
    return node;
}

BOOST_AUTO_TEST_CASE(ConstraintCMatrix)
{
    using namespace NuTo::Constraint;
    constexpr int numDofsTotal = 5;
    fakeit::Mock<NuTo::NodeBase> node = GetMockNode({3, 1});

    Constraints c;

    // 0th component of node --> dof = 3
    c.Add(eDofDisp, Equation({Term(node.get(), 0, 1)}));
    BOOST_CHECK(c.HaveChanged());
    BOOST_CHECK_EQUAL(c.GetNumEquations(eDofDisp), 1);

    auto m = c.BuildConstraintMatrix(eDofDisp, numDofsTotal).ConvertToFullMatrix();
    Eigen::MatrixXd expected(1, numDofsTotal);
    expected(0, 3) = 1;
    BoostUnitTest::CheckEigenMatrix(m, expected);

    // 1st component of node --> dof = 1
    c.Add(eDofDisp, Equation({Term(node.get(), 1, 1)}));
    BOOST_CHECK_EQUAL(c.GetNumEquations(eDofDisp), 2);
    expected.setZero(2, numDofsTotal);
    expected(0, 3) = 1;
    expected(1, 1) = 1;
    m = c.BuildConstraintMatrix(eDofDisp, numDofsTotal).ConvertToFullMatrix();
    BoostUnitTest::CheckEigenMatrix(m, expected);

    c.RemoveAll();
    BOOST_CHECK_EQUAL(c.GetNumEquations(eDofDisp), 0);
}

BOOST_AUTO_TEST_CASE(ConstraintRhs)
{
    using namespace NuTo::Constraint;
    Constraints c;
    c.Add(eDofDisp, Equation(1));
    c.Add(eDofDisp, Equation([](double time) { return time * 42; }));

    BoostUnitTest::CheckVector(c.GetRhs(eDofDisp, 0), std::vector<double>({1., 0.}), 2);
    BoostUnitTest::CheckVector(c.GetRhs(eDofDisp, 0.5), std::vector<double>({1., 21}), 2);
    BoostUnitTest::CheckVector(c.GetRhs(eDofDisp, 1), std::vector<double>({1., 42.}), 2);
}

BOOST_AUTO_TEST_CASE(ConstraintExchange)
{
    using namespace NuTo::Constraint;
    fakeit::Mock<NuTo::NodeBase> nodeA = GetMockNode({3});
    fakeit::Mock<NuTo::NodeBase> nodeB = GetMockNode({1});
    Constraints c;

    // 0th component of node --> dof = 3
    c.Add(eDofDisp, Equation({Term(nodeA.get(), 0, 1)}));
    auto m = c.BuildConstraintMatrix(eDofDisp, 5).ConvertToFullMatrix();
    Eigen::MatrixXd expected(1, 5);
    expected(0, 3) = 1;
    BoostUnitTest::CheckEigenMatrix(m, expected);

    c.ExchangeNodePtr(nodeA.get(), nodeB.get());
    expected.setZero();
    expected(0, 1) = 1;
    m = c.BuildConstraintMatrix(eDofDisp, 5).ConvertToFullMatrix();
    BoostUnitTest::CheckEigenMatrix(m, expected);
}
