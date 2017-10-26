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


// Mock<NuTo::NodeSimple> nodeMock;


NuTo::NodeSimple Node1(3);
NuTo::NodeSimple Node2(3);

// mock node, nodegroup, direction vector and lamdba function
class Helpers
{
public:
    Helpers()
        : node1(Node1)
        , node2(Node2)
    {
        nodes.Add(Node1);
        nodes.Add(Node2);
        someDirection << 1, 1, 1;
    }

    std::function<double(double)> rhsValueFunc(double value)
    {
        return [value](double) { return value; };
    }

    std::function<double(double)> sinFunction = [](double time) { return std::sin(time); };
    // std::function<double(double)> rhsOne = [](double) { return 1.; };
    Groups::Group<NodeSimple> nodes;
    NodeSimple& node1;
    NodeSimple& node2;
    Eigen::Vector3d someDirection;
    //    double reciprocalNorm = 1.0 / std::sqrt(3.0);
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
        expectedEquations.push_back(ConstraintPde::Equation(node, 0, rhsValueFunc(1.)));
    }
    BOOST_CHECK(eqs == expectedEquations);

    // with Group and lambda
    expectedEquations.clear();
    eqs = ConstraintPde::Component(nodes, {eDirection::X}, sinFunction);
    for (auto& node : nodes)
    {
        expectedEquations.push_back(ConstraintPde::Equation(node, 0, sinFunction));
    }
    BOOST_CHECK(eqs == expectedEquations);
}


BOOST_FIXTURE_TEST_CASE(Direction_test, Helpers)
{
    //    // with constant RHS
    //    auto eq = Constraint::Direction(node, someDirection);
    //    auto expectedEquation =
    //            Constraint::Equation({Constraint::Term(node, 0, reciprocalNorm), Constraint::Term(node, 1,
    //            reciprocalNorm),
    //                                  Constraint::Term(node, 2, reciprocalNorm)});
    //    BOOST_CHECK(eq == expectedEquation);

    //    // with lambda RHS
    //    eq = Constraint::Direction(node, someDirection, sinFunction);
    //    expectedEquation =
    //            Constraint::Equation({Constraint::Term(node, 0, reciprocalNorm), Constraint::Term(node, 1,
    //            reciprocalNorm),
    //                                  Constraint::Term(node, 2, reciprocalNorm)},
    //                                 sinFunction);
    //    BOOST_CHECK(eq == expectedEquation);

    //    // with Group
    //    std::vector<Constraint::Equation> expectedEquations;
    //    auto eqs = Constraint::Direction(nodes, someDirection);
    //    for (auto& node : nodes)
    //    {
    //        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, reciprocalNorm),
    //                                                          Constraint::Term(*node.second, 1, reciprocalNorm),
    //                                                          Constraint::Term(*node.second, 2, reciprocalNorm)}));
    //    }
    //    BOOST_CHECK(eqs == expectedEquations);

    //    // with Group and lambda
    //    expectedEquations.clear();
    //    eqs = Constraint::Direction(nodes, someDirection, sinFunction);
    //    for (auto& node : nodes)
    //    {
    //        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, reciprocalNorm),
    //                                                          Constraint::Term(*node.second, 1, reciprocalNorm),
    //                                                          Constraint::Term(*node.second, 2, reciprocalNorm)},
    //                                                         sinFunction));
    //    }
    //    BOOST_CHECK(eqs == expectedEquations);
}


// BOOST_FIXTURE_TEST_CASE(Value, Helpers)
//{
//    // with constant RHS
//    auto eq = Constraint::Value(node);
//    auto expectedEquation = Constraint::Equation({Constraint::Term(node, 0, 1.0)});
//    BOOST_CHECK(eq == expectedEquation);

//    // with lamdba RHS
//    eq = Constraint::Value(node, sinFunction);
//    expectedEquation = Constraint::Equation({Constraint::Term(node, 0, 1.0)}, sinFunction);
//    BOOST_CHECK(eq == expectedEquation);

//    // with Group
//    auto eqs = Constraint::Value(nodes);
//    std::vector<Constraint::Equation> expectedEquations;
//    for (auto& node : nodes)
//    {
//        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, 1.0)}));
//    }
//    BOOST_CHECK(eqs == expectedEquations);

//    // with Group
//    eqs = Constraint::Value(nodes, sinFunction);
//    expectedEquations.clear();
//    for (auto& node : nodes)
//    {
//        expectedEquations.push_back(Constraint::Equation({Constraint::Term(*node.second, 0, 1.0)}, sinFunction));
//    }
//    BOOST_CHECK(eqs == expectedEquations);
//}
