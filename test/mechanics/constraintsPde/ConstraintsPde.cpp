#include "BoostUnitTest.h"
#include "base/Exception.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/constraintsPde/Constraints.h"

BOOST_AUTO_TEST_CASE(ConstraintCMatrix)
{
    using namespace NuTo::ConstraintPde;
    constexpr int numDofsTotal = 5;
    
    NuTo::NodeSimple node(Eigen::Vector2d::Zero());
    node.SetDofNumber(0, 3);
    node.SetDofNumber(1, 1);

    Constraints c;

    NuTo::DofType dof("Erdös", 4);

    // 0th component of node --> dof = 3
    c.Add(dof, Equation(node, 0, [](double) {return 0;})); // Dof0 (3) * 1 = 0;
    BOOST_CHECK(c.HaveChanged());
    BOOST_CHECK_EQUAL(c.GetNumEquations(dof), 1);

    auto m = Eigen::MatrixXd(c.BuildConstraintMatrix(dof, numDofsTotal));
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(1, numDofsTotal);
    expected(0, 3) = 1;
    BoostUnitTest::CheckEigenMatrix(m, expected);

    // 1st component of node --> dof = 1
    c.Add(dof, Equation(node, 1, [](double) {return 42;})); // Dof1 (1) * 1 = 42;
    BOOST_CHECK_EQUAL(c.GetNumEquations(dof), 2);
    expected.setZero(2, numDofsTotal);
    expected(0, 3) = 1;
    expected(1, 1) = 1;
    m = Eigen::MatrixXd(c.BuildConstraintMatrix(dof, numDofsTotal));
    BoostUnitTest::CheckEigenMatrix(m, expected);
}

BOOST_AUTO_TEST_CASE(ConstraintRhs)
{
    using namespace NuTo::ConstraintPde;
    Constraints c;
    NuTo::NodeSimple dummyNode(Eigen::Vector2d::Zero());
    NuTo::DofType dof("Erdös", 4);
    c.Add(dof, Equation(dummyNode, 0, [](double){return 1;}));
    c.Add(dof, Equation(dummyNode, 1, [](double time) { return time * 42; }));

    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 0.0), Eigen::Vector2d(1, 0));
    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 0.5), Eigen::Vector2d(1, 21));
    BoostUnitTest::CheckEigenMatrix(c.GetRhs(dof, 1.0), Eigen::Vector2d(1, 42));
}

BOOST_AUTO_TEST_CASE(ConstraintEdgeCases)
{
    NuTo::NodeSimple dummyNode(Eigen::Vector2d::Zero()); // node with two components
    BOOST_CHECK_THROW(NuTo::ConstraintPde::Term(dummyNode, 42, 1), NuTo::Exception);
    BOOST_CHECK_THROW(NuTo::ConstraintPde::Term(dummyNode, 2, 1), NuTo::Exception);
    BOOST_CHECK_NO_THROW(NuTo::ConstraintPde::Term(dummyNode, 1, 1));
}
