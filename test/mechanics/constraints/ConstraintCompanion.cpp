#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include <cmath>
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/Group.h"

using namespace NuTo;
using namespace fakeit;

// Helper functions for comparisons
namespace NuTo
{
namespace Constraint
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
    auto ramp = Constraint::RhsRamp(42.0, 6174.0);
    BOOST_CHECK_EQUAL(ramp(0.0), 0.0);
    BOOST_CHECK_EQUAL(ramp(42.0), 6174.0);

    auto constant = Constraint::RhsConstant(42.0);
    BOOST_CHECK_EQUAL(constant(0.0), 42.0);
    BOOST_CHECK_EQUAL(constant(42.0), 42.0);
}


Mock<NuTo::NodeBase> nodeMock;

// mock node, nodegroup, direction vector and lamdba function
class Helpers
{
public:
    Helpers()
        : node(nodeMock.get())
    {
        When(Method(nodeMock, GetNumDofs)).AlwaysReturn(3);
        nodes.AddMember(0, &nodeMock.get());
        nodes.AddMember(1, &nodeMock.get());
        someDirection << 1, 1, 1;
    }
    std::function<double(double)> sinFunction = [](double time) { return std::sin(time); };
    Group<NodeBase> nodes;
    NodeBase& node;
    Eigen::Vector3d someDirection;
    double reciprocalNorm = 1.0 / std::sqrt(3.0);
};


BOOST_FIXTURE_TEST_CASE(Component, Helpers)
{
    // with constant RHS
    auto eqs = Constraint::Component(node, {eDirection::X});
    auto expectedEquation = Constraint::Equation({Constraint::Term(node, 0, 1.0)});
    BOOST_CHECK(eqs[0] == expectedEquation);

    // with lambda RHS
    eqs = Constraint::Component(node, {eDirection::X}, sinFunction);
    expectedEquation = Constraint::Equation({Constraint::Term(node, 0, 1.0)}, sinFunction);
    BOOST_CHECK(eqs[0] == expectedEquation);

    // with Group
    std::vector<Constraint::Equation> expectedEquations;
    eqs = Constraint::Component(nodes, {eDirection::X});
    for (auto& node : nodes)
    {
        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, 1.0)}));
    }
    BOOST_CHECK(eqs == expectedEquations);

    // with Group and lambda
    expectedEquations.clear();
    eqs = Constraint::Component(nodes, {eDirection::X}, sinFunction);
    for (auto& node : nodes)
    {
        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, 1.0)}, sinFunction));
    }
    BOOST_CHECK(eqs == expectedEquations);
}


BOOST_FIXTURE_TEST_CASE(Direction, Helpers)
{
    // with constant RHS
    auto eq = Constraint::Direction(node, someDirection);
    auto expectedEquation =
            Constraint::Equation({Constraint::Term(node, 0, reciprocalNorm), Constraint::Term(node, 1, reciprocalNorm),
                                  Constraint::Term(node, 2, reciprocalNorm)});
    BOOST_CHECK(eq == expectedEquation);

    // with lambda RHS
    eq = Constraint::Direction(node, someDirection, sinFunction);
    expectedEquation =
            Constraint::Equation({Constraint::Term(node, 0, reciprocalNorm), Constraint::Term(node, 1, reciprocalNorm),
                                  Constraint::Term(node, 2, reciprocalNorm)},
                                 sinFunction);
    BOOST_CHECK(eq == expectedEquation);

    // with Group
    std::vector<Constraint::Equation> expectedEquations;
    auto eqs = Constraint::Direction(nodes, someDirection);
    for (auto& node : nodes)
    {
        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, reciprocalNorm),
                                                          Constraint::Term(*node.second, 1, reciprocalNorm),
                                                          Constraint::Term(*node.second, 2, reciprocalNorm)}));
    }
    BOOST_CHECK(eqs == expectedEquations);

    // with Group and lambda
    expectedEquations.clear();
    eqs = Constraint::Direction(nodes, someDirection, sinFunction);
    for (auto& node : nodes)
    {
        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, reciprocalNorm),
                                                          Constraint::Term(*node.second, 1, reciprocalNorm),
                                                          Constraint::Term(*node.second, 2, reciprocalNorm)},
                                                         sinFunction));
    }
    BOOST_CHECK(eqs == expectedEquations);
}


BOOST_FIXTURE_TEST_CASE(Value, Helpers)
{
    // with constant RHS
    auto eq = Constraint::Value(node);
    auto expectedEquation = Constraint::Equation({Constraint::Term(node, 0, 1.0)});
    BOOST_CHECK(eq == expectedEquation);

    // with lamdba RHS
    eq = Constraint::Value(node, sinFunction);
    expectedEquation = Constraint::Equation({Constraint::Term(node, 0, 1.0)}, sinFunction);
    BOOST_CHECK(eq == expectedEquation);

    // with Group
    auto eqs = Constraint::Value(nodes);
    std::vector<Constraint::Equation> expectedEquations;
    for (auto& node : nodes)
    {
        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, 1.0)}));
    }
    BOOST_CHECK(eqs == expectedEquations);

    // with Group
    eqs = Constraint::Value(nodes, sinFunction);
    expectedEquations.clear();
    for (auto& node : nodes)
    {
        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, 1.0)}, sinFunction));
    }
    BOOST_CHECK(eqs == expectedEquations);
}
