#include "BoostUnitTest.h"
#include "mechanics/constraintsPde/ConstraintCompanion.h"

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
    if (std::abs(a.GetCoefficient() - b.GetCoefficient()) > 1.e-6)
        return false;
    if (a.GetComponent() != b.GetComponent())
        return false;
    return true;
}

bool operator==(const Equation& a, const Equation& b)
{
    if (std::abs(a.GetRhs(0.0) - b.GetRhs(0.0)) > 1.e-6)
        return false;
    if (std::abs(a.GetRhs(42.0) - b.GetRhs(42.0)) > 1.e-6)
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
    someDirection << 2, -3, 4;

    // single equation with constant RHS
    auto eq = ConstraintPde::Direction(node1, someDirection);
    auto expectedEquation = ConstraintPde::Equation(node1, 2, rhsValueFunc(0.0));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 0, -0.5));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 1, 0.75));

    BOOST_CHECK(eq == expectedEquation);

    // single equation with lambda RHS
    eq = ConstraintPde::Direction(node1, someDirection, sinFunction);
    expectedEquation = ConstraintPde::Equation(node1, 2, rhsValueAdjustedFunc(sinFunction, std::sqrt(29.) / 4.));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 0, -0.5));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 1, 0.75));

    BOOST_CHECK(eq == expectedEquation);

    // group with constant RHS
    std::vector<ConstraintPde::Equation> expectedEquations;
    auto eqs = ConstraintPde::Direction(nodes, someDirection);
    for (auto& node : nodes)
    {
        expectedEquation = ConstraintPde::Equation(node, 2, rhsValueFunc(0.0));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 0, -0.5));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 1, 0.75));
        expectedEquations.push_back(expectedEquation);
    }
    BOOST_CHECK(eqs == expectedEquations);

    // group with lambda RHS
    expectedEquations.clear();
    eqs = ConstraintPde::Direction(nodes, someDirection, sinFunction);
    for (auto& node : nodes)
    {
        expectedEquation = ConstraintPde::Equation(node, 2, rhsValueAdjustedFunc(sinFunction, std::sqrt(29.) / 4.));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 0, -0.5));
        expectedEquation.AddTerm(ConstraintPde::Term(node, 1, 0.75));
        expectedEquations.push_back(expectedEquation);
    }
    BOOST_CHECK(eqs == expectedEquations);


    // Tests if terms with a coefficient of 0 are ignored correctly

    someDirection << 0., 3., -4.;
    eq = ConstraintPde::Direction(node1, someDirection, sinFunction);
    expectedEquation = ConstraintPde::Equation(node1, 2, rhsValueAdjustedFunc(sinFunction, -5. / 4.));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 1, 0.75));

    BOOST_CHECK(eq == expectedEquation);
    BOOST_CHECK_EQUAL(eq.GetTerms().size(), 2);


    someDirection << 3., 0., -4.;
    eq = ConstraintPde::Direction(node1, someDirection, sinFunction);
    expectedEquation = ConstraintPde::Equation(node1, 2, rhsValueAdjustedFunc(sinFunction, -5. / 4.));
    expectedEquation.AddTerm(ConstraintPde::Term(node1, 0, 0.75));

    BOOST_CHECK(eq == expectedEquation);
    BOOST_CHECK_EQUAL(eq.GetTerms().size(), 2);


    someDirection << 0., 0., 0.;
    BOOST_CHECK_THROW(ConstraintPde::Direction(node1, someDirection, sinFunction), Exception);
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
