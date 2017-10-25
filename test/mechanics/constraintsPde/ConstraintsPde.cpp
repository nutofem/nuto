#include "BoostUnitTest.h"
#include "base/Exception.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/constraintsPde/Constraints.h"

using namespace NuTo;
using namespace ConstraintPde;

const DofType dof("Erdös", 4);
auto rhs = [](double) { return 42; };

BOOST_AUTO_TEST_CASE(ConstraintCMatrix)
{
    constexpr int numDofsTotal = 5;

    NodeSimple node(Eigen::Vector2d::Zero());
    node.SetDofNumber(0, 3);
    node.SetDofNumber(1, 1);

    Constraints c;

    // 0th component of node --> dof = 3
    c.Add(dof, Equation(node, 0, [](double) { return 0; })); // Dof0 (3) * 1 = 0;
    BOOST_CHECK(c.HaveChanged());
    BOOST_CHECK_EQUAL(c.GetNumEquations(dof), 1);

    auto m = Eigen::MatrixXd(c.BuildConstraintMatrix(dof, numDofsTotal));
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(1, numDofsTotal);
    expected(0, 3) = 1;
    BoostUnitTest::CheckEigenMatrix(m, expected);

    // 1st component of node --> dof = 1
    c.Add(dof, Equation(node, 1, [](double) { return 42; })); // Dof1 (1) * 1 = 42;
    BOOST_CHECK_EQUAL(c.GetNumEquations(dof), 2);
    expected.setZero(2, numDofsTotal);
    expected(0, 3) = 1;
    expected(1, 1) = 1;
    m = Eigen::MatrixXd(c.BuildConstraintMatrix(dof, numDofsTotal));
    BoostUnitTest::CheckEigenMatrix(m, expected);
}

BOOST_AUTO_TEST_CASE(ConstraintRhs)
{
    Constraints c;
    NodeSimple dummyNode(Eigen::Vector2d::Zero());
    c.Add(dof, Equation(dummyNode, 0, [](double) { return 1; }));
    c.Add(dof, Equation(dummyNode, 1, [](double time) { return time * 42; }));

    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 0.0), Eigen::Vector2d(1, 0));
    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 0.5), Eigen::Vector2d(1, 21));
    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 1.0), Eigen::Vector2d(1, 42));
}

BOOST_AUTO_TEST_CASE(ConstraintUnavailableComponent)
{
    NodeSimple dummyNode(Eigen::Vector2d::Zero()); // node with two components
    BOOST_CHECK_THROW(Term(dummyNode, 42, 1), Exception);
    BOOST_CHECK_THROW(Term(dummyNode, 2, 1), Exception);
    BOOST_CHECK_NO_THROW(Term(dummyNode, 1, 1));
}

BOOST_AUTO_TEST_CASE(ConstraintDoubleConstrained)
{
    Constraints constraints;
    NodeSimple node(Eigen::Vector2d::Zero());

    constraints.Add(dof, {node, 0, rhs});
    BOOST_CHECK_THROW(constraints.Add(dof, {node, 0, rhs}), Exception);
}

BOOST_AUTO_TEST_CASE(ConstraintTwoDependentDofsInOneEquation)
{
    ConstraintPde::Constraints constraints;
    NodeSimple node(Eigen::Vector2d::Zero());

    constraints.Add(dof, {node, 0, rhs}); // node.dof0 * 1.0 = rhs

    ConstraintPde::Equation equation(node, 1, rhs);
    equation.AddTerm({node, 0, 0.4}); // node.dof1 * 1.0 + node.dof0 * 0.4 = rhs
    BOOST_CHECK_THROW(constraints.Add(dof, equation), Exception);
}
