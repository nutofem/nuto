#include <iostream>
#include "BoostUnitTest.h"
#include "nuto/base/Exception.h"
#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/constraints/Constraints.h"

using namespace NuTo;
using namespace Constraint;

const DofType dof("Erdös", 4);
auto rhs = [](double) { return 42; };

BOOST_AUTO_TEST_CASE(ConstraintUnnumbered)
{
    DofNode node(Eigen::Vector2d::Zero());
    Constraints c;
    c.Add(dof, Equation(node, 0, rhs));
    // the dofs are not numbered.
    BOOST_CHECK_THROW(c.BuildUnitConstraintMatrix(dof, 2), Exception);
}

BOOST_AUTO_TEST_CASE(ConstraintCMatrix)
{
    DofNode node0(0);
    DofNode node1(0);
    DofNode node2(0);
    DofNode node3(0);

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

    BOOST_CHECK_NO_THROW(c.BuildUnitConstraintMatrix(dof, 4));
}

BOOST_AUTO_TEST_CASE(ConstraintCMatrixInteracting)
{
    DofNode node0(0);
    DofNode node1(0);
    DofNode node2(0);
    DofNode node3(0);
    DofNode node4(0);

    /*
     *     n0 ---- n1 ---- n2 ---- n3 --- n4 --- n5
     *         n1(dep)           = rhs;
     *    with n4(dep) + 42 * n2 = rhs;
     *
     *                (n0)
     *                (n2)
     *  [ 0  0  0  0] (n3) + [ 1  0 ] (n1) = rhs;
     *  [ 0 42  0  0] (n5) + [ 0  1 ] (n4) = rhs;
     *
     */
    Eigen::MatrixXd cmatUnitExpected = Eigen::MatrixXd::Zero(6, 4);
    cmatUnitExpected(0, 0) = 1;
    cmatUnitExpected(2, 1) = 1;
    cmatUnitExpected(3, 2) = 1;
    cmatUnitExpected(4, 1) = -42;
    cmatUnitExpected(5, 3) = 1;

    node0.SetDofNumber(0, 0);
    node1.SetDofNumber(0, 1);
    node2.SetDofNumber(0, 2);
    node3.SetDofNumber(0, 3);
    node4.SetDofNumber(0, 4);
    node5.SetDofNumber(0, 5);

    Equation noninteractingEquation(node1, 0, rhs);

    Equation interactingEquation(node4, 0, rhs);
    interactingEquation.AddIndependentTerm({node2, 0, 42});

    Constraints c;
    c.Add(dof, noninteractingEquation);
    c.Add(dof, interactingEquation);

    Eigen::MatrixXd cmatUnit = c.BuildUnitConstraintMatrix(dof, 6);

    std::cout << "cmatUnitExpected\n" << cmatUnitExpected << std::endl;
    std::cout << "cmatUnit\n" << cmatUnit << std::endl;
    BoostUnitTest::CheckEigenMatrix(cmatUnit, cmatUnitExpected);
}


BOOST_AUTO_TEST_CASE(ConstraintRhs)
{
    Constraints c;
    DofNode dummyNode(Eigen::Vector2d::Zero());
    c.Add(dof, Equation(dummyNode, 0, [](double) { return 1; }));
    c.Add(dof, Equation(dummyNode, 1, [](double time) { return time * 42; }));

    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 0.0), Eigen::Vector2d(1, 0));
    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 0.5), Eigen::Vector2d(1, 21));
    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 1.0), Eigen::Vector2d(1, 42));
}

BOOST_AUTO_TEST_CASE(ConstraintUnavailableComponent)
{
    DofNode dummyNode(Eigen::Vector2d::Zero()); // node with two components
    BOOST_CHECK_THROW(Term(dummyNode, 42, 1), Exception);
    BOOST_CHECK_THROW(Term(dummyNode, 2, 1), Exception);
    BOOST_CHECK_NO_THROW(Term(dummyNode, 1, 1));
}

BOOST_AUTO_TEST_CASE(ConstraintDoubleConstrained)
{
    Constraints constraints;
    DofNode node(Eigen::Vector2d::Zero());

    constraints.Add(dof, {node, 0, rhs});
    BOOST_CHECK_THROW(constraints.Add(dof, {node, 0, rhs}), Exception);
}

BOOST_AUTO_TEST_CASE(ConstraintTwoDependentDofsInOneEquation)
{
    Constraint::Constraints constraints;
    DofNode node(Eigen::Vector2d::Zero());

    constraints.Add(dof, {node, 0, rhs}); // node.dof0 * 1.0 = rhs

    Constraint::Equation equation(node, 1, rhs);
    equation.AddIndependentTerm({node, 0, 0.4}); // node.dof1 * 1.0 + node.dof0 * 0.4 = rhs
    BOOST_CHECK_THROW(constraints.Add(dof, equation), Exception);
}

BOOST_AUTO_TEST_CASE(ConstraintEmptyEquations)
{
    Constraint::Constraints constraints;
    // passing an empty vector is _most likely_ not intended and should throw
    BOOST_CHECK_THROW(constraints.Add(dof, {}), Exception);
}
