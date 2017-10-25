#include "BoostUnitTest.h"
#include "base/Exception.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/constraintsPde/Constraints.h"

using namespace NuTo;
using namespace ConstraintPde;

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
