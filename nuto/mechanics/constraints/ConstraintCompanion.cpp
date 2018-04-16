#include "nuto/base/Exception.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/constraints/Equation.h"


namespace NuTo
{
namespace Constraint
{

RhsFunction RhsRamp(double timeEnd, double valueEnd)
{
    if (timeEnd == 0)
        throw Exception(__PRETTY_FUNCTION__, "The first parameter must be != 0. If you wanted a constant rhs, use "
                                             "Constraint::RhsConstant(constant).");
    return [=](double time) { return valueEnd * time / timeEnd; };
}


std::vector<Equation> Component(const DofNode& node, std::vector<eDirection> directions, double value)
{
    return Component(node, directions, RhsConstant(value));
}


std::vector<Equation> Component(const DofNode& node, std::vector<eDirection> directions, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto direction : directions)
    {
        int component = ToComponentIndex(direction);
        eqs.emplace_back(node, component, rhs);
    }
    return eqs;
}


std::vector<Equation> Component(const Group<DofNode>& nodes, std::vector<eDirection> directions, double value)
{
    return Component(nodes, directions, RhsConstant(value));
}


std::vector<Equation> Component(const Group<DofNode>& nodes, std::vector<eDirection> directions, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (const auto& node : nodes)
    {
        auto tmpEqs = Component(node, directions, rhs);
        eqs.insert(eqs.end(), tmpEqs.begin(), tmpEqs.end());
    }
    return eqs;
}


Equation Direction(const DofNode& node, Eigen::VectorXd direction, RhsFunction rhs)
{

    int maxComponentIndex = -1;
    if (direction.cwiseAbs().maxCoeff(&maxComponentIndex) < 1.e-6)
        throw Exception(__PRETTY_FUNCTION__,
                        "Your direction vector is composed of zeros only! The direction is unspecified!");

    // Normalization necessary, otherwise rhs depends on length of the direction vector
    direction.normalize();

    // first nonzero direction component defines the dependent dof
    // Lambda corrects original rhs funtion
    Equation e(node, maxComponentIndex,
               [=](double time) -> double { return rhs(time) / direction[maxComponentIndex]; });

    // add terms for all non zero direction components
    for (int iComponent = 0; iComponent < direction.rows(); ++iComponent)
        if (std::abs(direction[iComponent]) > 0 && iComponent != maxComponentIndex)
            e.AddTerm(Term(node, iComponent, -direction[iComponent] / direction[maxComponentIndex]));


    return e;
}

Equation Direction(const DofNode& node, Eigen::VectorXd direction, double value)
{
    return Direction(node, direction, RhsConstant(value));
}

std::vector<Equation> Direction(const Group<DofNode>& nodes, Eigen::VectorXd direction, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto& node : nodes)
        eqs.push_back(Direction(node, direction, rhs));
    return eqs;
}

std::vector<Equation> Direction(const Group<DofNode>& nodes, Eigen::VectorXd direction, double value)
{
    return Direction(nodes, direction, RhsConstant(value));
}

Equation Value(const DofNode& node, double value)
{
    return Value(node, RhsConstant(value));
}

Equation Value(const DofNode& node, RhsFunction rhs)
{
    if (node.GetNumValues() != 1)
        throw Exception(__PRETTY_FUNCTION__, "This function is meant to be used with single value nodes only");
    return Component(node, {eDirection::X}, rhs)[0];
}

std::vector<Equation> Value(const Group<DofNode>& nodes, double value)
{
    std::vector<Equation> eqs;
    for (auto& node : nodes)
        eqs.push_back(Value(node, value));
    return eqs;
}

std::vector<Equation> Value(const Group<DofNode>& nodes, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto& node : nodes)
        eqs.push_back(Value(node, rhs));
    return eqs;
}

} /* ConstraintPde */
} /* NuTo */
