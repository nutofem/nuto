#include "BoostUnitTest.h"
#include <fakeit.hpp>
//#include <cmath>
//#include "mechanics/groups/Group.h"
#include "mechanics/constraintsPde/ConstraintCompanion.h"

using namespace fakeit;
using namespace NuTo;
using namespace NuTo::ConstraintPde;


//// Helper functions for comparisons
namespace NuTo
{
namespace ConstraintPde
{
bool operator==(const Term& a, const Term& b)
{
    if (&(a.GetNode()) != &(b.GetNode()))
        return false;
    if (a.GetCoefficient() != b.GetCoefficient())
        return false;
    if (a.GetComponent() != b.GetComponent())
        return false;
    return true;
}

bool operator==(const Equation& a, const Equation& b)
{
    if (a.GetRhs(0.0) != b.GetRhs(0.0))
        return false;
    if (a.GetRhs(42.0) != b.GetRhs(42.0))
        return false;
    if (a.GetTerms() != b.GetTerms())
        return false;
    return true;
}

} // Constraint
} // NuTo


// quite the stupid tests, but oh well...
BOOST_AUTO_TEST_CASE(PredefinedFunctions)
{
    auto ramp = ConstraintPde::RhsRamp(42.0, 6174.0);
    BOOST_CHECK_EQUAL(ramp(0.0), 0.0);
    BOOST_CHECK_EQUAL(ramp(42.0), 6174.0);

    auto constant = ConstraintPde::RhsConstant(42.0);
    BOOST_CHECK_EQUAL(constant(0.0), 42.0);
    BOOST_CHECK_EQUAL(constant(42.0), 42.0);
}


// mock node, nodegroup, direction vector and lamdba function
class Helpers
{
public:
    Helpers()
    {
        nodes.Add(node1);
        nodes.Add(node2);
    }

    std::function<double(double)> rhsValueFunc(double value)
    {
        return [value](double) { return value; };
    }
    std::function<double(double)> rhsValueAdjustedFunc(std::function<double(double)> f, double value)
    {
        return [f, value](double time) { return f(time) * value; };
    }

    std::function<double(double)> sinFunction = [](double time) { return std::sin(time); };
    Groups::Group<NodeSimple> nodes;
    NodeSimple node1 = NodeSimple(Eigen::VectorXd::Ones(3));
    NodeSimple node2 = NodeSimple(Eigen::VectorXd::Ones(3));
    double reciprocalNorm = 1.0 / std::sqrt(3.0);
};


BOOST_FIXTURE_TEST_CASE(Component_test, Helpers)
{
    // with constant RHS
    auto eqs = ConstraintPde::Component(node1, {eDirection::X}, 1.0);
    auto expectedEquation = ConstraintPde::Equation(node1, 0, rhsValueFunc(1.));
    BOOST_CHECK(eqs[0] == expectedEquation);

    // with lambda RHS
    eqs = ConstraintPde::Component(node1, {eDirection::X}, sinFunction);
    expectedEquation = ConstraintPde::Equation(node1, 0, sinFunction);
    BOOST_CHECK(eqs[0] == expectedEquation);

    // with Group
    std::vector<ConstraintPde::Equation> expectedEquations;
    eqs = ConstraintPde::Component(nodes, {eDirection::X}, 1.);
    for (auto& node : nodes)
    {
        expectedEquations.emplace_back(node, 0, rhsValueFunc(1.));
    }
    BOOST_CHECK(eqs == expectedEquations);

    // with Group and lambda
    expectedEquations.clear();
    eqs = ConstraintPde::Component(nodes, {eDirection::X}, sinFunction);
    for (auto& node : nodes)
    {
        expectedEquations.emplace_back(node, 0, sinFunction);
    }
    BOOST_CHECK(eqs == expectedEquations);
}


BOOST_FIXTURE_TEST_CASE(Direction_test, Helpers)
{
    Eigen::Vector3d someDirection;
    someDirection << 1, 1, 1;

    // single equation with constant RHS
    auto eq = ConstraintPde::Direction(node1, someDirection);
    auto expectedEquation = ConstraintPde::Equation(node1, 0, rhsValueFunc(0.0));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 1, -1.));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 2, -1.));

    BOOST_CHECK(eq == expectedEquation);

    // single equation with lambda RHS
    eq = ConstraintPde::Direction(node1, someDirection, sinFunction);
    expectedEquation = ConstraintPde::Equation(node1, 0, rhsValueAdjustedFunc(sinFunction, 1. / reciprocalNorm));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 1, -1.));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 2, -1.));

    BOOST_CHECK(eq == expectedEquation);

    // group with constant RHS
    std::vector<ConstraintPde::Equation> expectedEquations;
    auto eqs = ConstraintPde::Direction(nodes, someDirection);
    for (auto& node : nodes)
    {
        expectedEquation = ConstraintPde::Equation(node, 0, rhsValueFunc(0.0));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 1, -1.));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 2, -1.));
        expectedEquations.push_back(expectedEquation);
    }
    BOOST_CHECK(eqs == expectedEquations);

    // group with lambda RHS
    expectedEquations.clear();
    eqs = ConstraintPde::Direction(nodes, someDirection, sinFunction);
    for (auto& node : nodes)
    {
        expectedEquation = ConstraintPde::Equation(node, 0, rhsValueAdjustedFunc(sinFunction, 1. / reciprocalNorm));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 1, -1.));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 2, -1.));
        expectedEquations.push_back(expectedEquation);
    }
    BOOST_CHECK(eqs == expectedEquations);
}


BOOST_FIXTURE_TEST_CASE(Value_test, Helpers)
{

    BOOST_CHECK_THROW(ConstraintPde::Value(node1, 42.), Exception);
    node1 = NodeSimple(Eigen::VectorXd::Ones(1));
    node2 = NodeSimple(Eigen::VectorXd::Ones(1));


    // with constant RHS
    auto eq = ConstraintPde::Value(node1, 42.);
    auto expectedEquation = ConstraintPde::Equation(node1, 0, rhsValueFunc(42.0));
    BOOST_CHECK(eq == expectedEquation);

    // with lamdba RHS
    eq = ConstraintPde::Value(node1, sinFunction);
    expectedEquation = ConstraintPde::Equation(node1, 0, sinFunction);
    BOOST_CHECK(eq == expectedEquation);

    // with Group
    auto eqs = ConstraintPde::Value(nodes, 42.0);
    std::vector<ConstraintPde::Equation> expectedEquations;
    for (auto& node : nodes)
    {
        expectedEquations.emplace_back(node, 0, rhsValueFunc(42.0));
    }
    BOOST_CHECK(eqs == expectedEquations);

    // with Group
    eqs = ConstraintPde::Value(nodes, sinFunction);
    expectedEquations.clear();
    for (auto& node : nodes)
    {
        expectedEquations.emplace_back(node, 0, sinFunction);
    }
    BOOST_CHECK(eqs == expectedEquations);
}
