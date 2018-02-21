#include "BoostUnitTest.h"
#include "base/Exception.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/constraints/Constraints.h"

using namespace NuTo;
using namespace Constraint;

const DofType dof("Erdös", 4);
auto rhs = [](double) { return 42; };

BOOST_AUTO_TEST_CASE(ConstraintUnnumbered)
{
    NodeSimple node(Eigen::Vector2d::Zero());
    Constraints c;
    c.Add(dof, Equation(node, 0, rhs));
    // the dofs are not numbered.
    BOOST_CHECK_THROW(c.BuildConstraintMatrix(dof, 1), Exception);
}

BOOST_AUTO_TEST_CASE(ConstraintCMatrix)
{
    NodeSimple node0(0);
    NodeSimple node1(0);
    NodeSimple node2(0);
    NodeSimple node3(0);

    /*
     *  n0 ---- n1 ---- n2 ---- n3
     *  //                      //
     * (fix, eq0)          (also fix, eq1)
     */

    Constraints c;
    c.Add(dof, Equation(node0, 0, rhs));
    c.Add(dof, Equation(node3, 0, rhs));

    // unconstrained
    node1.SetDofNumber(0, 0);
    node2.SetDofNumber(0, 1);

    // constrained, order 2,3 matters.
    node0.SetDofNumber(0, 2);
    node3.SetDofNumber(0, 3);

    BOOST_CHECK_NO_THROW(c.BuildConstraintMatrix(dof, 2));

    // provide invalid numbering of dependent dofs
    node0.SetDofNumber(0, 3);
    node3.SetDofNumber(0, 2);

    // this will not build an identity matrix at the end of CMat but
    // 0  1
    // 1  0
    // which is wrong.
    BOOST_CHECK_THROW(c.BuildConstraintMatrix(dof, 2), Exception);
}

BOOST_AUTO_TEST_CASE(ConstraintCMatrixInteracting)
{
    NodeSimple node0(0);
    NodeSimple node1(0);
    NodeSimple node2(0);
    NodeSimple node3(0);
    NodeSimple node4(0);

    /*
     *     n0 ---- n1 ---- n2 ---- n3 --- n4
     *         n3(dep)           = rhs;
     *    with n4(dep) + 42 * n0 = rhs;
     *
     *              (n0)
     *  [ 0  0  0 ] (n1) + [ 1  0 ] (n3) = rhs;
     *  [42  0  0 ] (n2) + [ 0  1 ] (n4) = rhs;
     *
     */
    Eigen::MatrixXd cmatExpected = Eigen::MatrixXd::Zero(2, 3);
    cmatExpected(1, 0) = 42;

    node0.SetDofNumber(0, 0);
    node1.SetDofNumber(0, 1);
    node2.SetDofNumber(0, 2);
    node3.SetDofNumber(0, 3);
    node4.SetDofNumber(0, 4);

    Equation noninteractingEquation(node3, 0, rhs);

    Equation interactingEquation(node4, 0, rhs);
    interactingEquation.AddTerm({node0, 0, 42});

    Constraints c;
    c.Add(dof, noninteractingEquation);
    c.Add(dof, interactingEquation);

    Eigen::MatrixXd cmat = c.BuildConstraintMatrix(dof, 3);

    BoostUnitTest::CheckEigenMatrix(cmat, cmatExpected);
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
    Constraint::Constraints constraints;
    NodeSimple node(Eigen::Vector2d::Zero());

    constraints.Add(dof, {node, 0, rhs}); // node.dof0 * 1.0 = rhs

    Constraint::Equation equation(node, 1, rhs);
    equation.AddTerm({node, 0, 0.4}); // node.dof1 * 1.0 + node.dof0 * 0.4 = rhs
    BOOST_CHECK_THROW(constraints.Add(dof, equation), Exception);
}

BOOST_AUTO_TEST_CASE(ConstraintEmptyEquations)
{
    Constraint::Constraints constraints;
    // passing an empty vector is _most likely_ not intended and should throw
    BOOST_CHECK_THROW(constraints.Add(dof, {}), Exception);
}
