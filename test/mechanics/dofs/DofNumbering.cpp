#include "BoostUnitTest.h"

#include "mechanics/dofs/DofNumbering.h"

// wrap algorithm to provide range-like interface
bool IsPermutation(std::vector<int> v0, std::vector<int> v1)
{
    return std::is_permutation(v0.begin(), v0.end(), v1.begin(), v1.end());
}

using namespace NuTo;

auto rhs = [](double) { return 42; };
const DofType d("stuff", 2);

BOOST_AUTO_TEST_CASE(DofNumberingTest)
{
    // Add three 2d nodes, one of them constrained in x and y.
    NodeSimple nodeUnconstrained0(Eigen::Vector2d::Zero());
    NodeSimple nodeUnconstrained1(Eigen::Vector2d::Zero());
    NodeSimple nodeConstrained(Eigen::Vector2d::Zero());

    Constraint::Constraints constraints;

    Groups::Group<NodeSimple> group({nodeConstrained, nodeUnconstrained0, nodeUnconstrained1});

    constraints.Add(d, {nodeConstrained, 1, rhs}); // 1st equation for component 1 --> expected number 4
    constraints.Add(d, {nodeConstrained, 0, rhs}); // 2nd equation for component 0 --> expected number 5

    auto dofInfo = DofNumbering::Build(group, d, constraints);
    BOOST_CHECK_EQUAL(dofInfo.numDependentDofs[d], 2);
    BOOST_CHECK_EQUAL(dofInfo.numIndependentDofs[d], 4);

    // Check if unconstrained nodes have dof number [0..3], ordering does not matter
    std::vector<int> dofNumbersUnconstrained{nodeUnconstrained0.GetDofNumber(0), nodeUnconstrained0.GetDofNumber(1),
                                             nodeUnconstrained1.GetDofNumber(0), nodeUnconstrained1.GetDofNumber(1)};
    BOOST_CHECK(IsPermutation(dofNumbersUnconstrained, {0, 1, 2, 3}));

    // Check if constrained nodes have dof number [5, 4], ordering has to be in the order of the equation definition
    BOOST_CHECK_EQUAL(nodeConstrained.GetDofNumber(0), 5);
    BOOST_CHECK_EQUAL(nodeConstrained.GetDofNumber(1), 4);
}
