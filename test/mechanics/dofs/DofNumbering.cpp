#include "BoostUnitTest.h"

#include "mechanics/dofs/DofNumbering.h"
#include <set>

using namespace NuTo;

auto rhs = [](double) { return 42; };

const DofType d("stuff", 2);

BOOST_AUTO_TEST_CASE(DofNumberingTest)
{
    /* Add three 2d nodes, one of them constrained in x and y.
     * Check if unconstrained nodes have dof number [0..3] and
     * constrained nodes have dof number [4..5].
     */
    NodeSimple nodeUnconstrained0(Eigen::Vector2d::Zero());
    NodeSimple nodeUnconstrained1(Eigen::Vector2d::Zero());
    NodeSimple nodeConstrained(Eigen::Vector2d::Zero());

    ConstraintPde::Constraints constraints;

    Groups::Group<NodeSimple> group({nodeConstrained, nodeUnconstrained0, nodeUnconstrained1});

    constraints.Add(d, {nodeConstrained, 1, rhs});
    constraints.Add(d, {nodeConstrained, 0, rhs});

    auto dofInfo = DofNumbering::Build(group, d, constraints);
    BOOST_CHECK_EQUAL(dofInfo.numDependentDofs[d], 2);
    BOOST_CHECK_EQUAL(dofInfo.numIndependentDofs[d], 4);

    std::set<int> dofNumbers{nodeUnconstrained0.GetDofNumber(0), nodeUnconstrained0.GetDofNumber(1),
                             nodeUnconstrained1.GetDofNumber(0), nodeUnconstrained1.GetDofNumber(1)};
    BOOST_CHECK_EQUAL(dofNumbers.size(), 4);
    BOOST_CHECK_EQUAL(*dofNumbers.rbegin(), 3);

    dofNumbers.insert(nodeConstrained.GetDofNumber(0));
    dofNumbers.insert(nodeConstrained.GetDofNumber(1));
    BOOST_CHECK_EQUAL(dofNumbers.size(), 6);
    BOOST_CHECK_EQUAL(*dofNumbers.rbegin(), 5);

    auto cmat = constraints.BuildConstraintMatrix(d, 6);
    Eigen::MatrixXd identityBlock = cmat.block(0, 4, 2, 2);
    BoostUnitTest::CheckEigenMatrix(identityBlock, Eigen::Matrix2d::Identity());
}


