#include "mechanics/constraintsPde/ConstraintCompanion.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/constraintsPde/Equation.h"

namespace NuTo
{
namespace ConstraintPde
{

std::vector<Equation> Component(const NodeSimple& node, std::vector<eDirection> directions, double value)
{
    return Component(node, directions, RhsConstant(value));
}


std::vector<Equation> Component(const NodeSimple& node, std::vector<eDirection> directions, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto direction : directions)
    {
        int component = ToComponentIndex(direction);
        eqs.emplace_back(node, component, rhs);
    }
    return eqs;
}


std::vector<Equation> Component(const Groups::Group<NodeSimple>& nodes, std::vector<eDirection> directions,
                                double value)
{
    return Component(nodes, directions, RhsConstant(value));
}


std::vector<Equation> Component(const Groups::Group<NodeSimple>& nodes, std::vector<eDirection> directions,
                                RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (const auto& node : nodes)
    {
        auto tmpEqs = Component(node, directions, rhs);
        eqs.insert(eqs.end(), tmpEqs.begin(), tmpEqs.end());
    }
    return eqs;
}


Equation Direction(const NodeSimple& node, Eigen::VectorXd direction, RhsFunction rhs)
{

    assert(direction != Eigen::VectorXd::Zero(direction.rows()) && "the direction vector contains only zeros");

    // Normalization necessary, otherwise rhs depends on length of the direction vector
    direction.normalize();

    int firstNonZeroComp = -1;
    for (int iComponent = 0; iComponent < direction.rows(); ++iComponent)
    {
        if (std::abs(direction[iComponent]) > 0)
        {
            firstNonZeroComp = iComponent;
            break;
        }
    }

    // first nonzero direction component defines the dependent dof
    // Lambda corrects original rhs funtion
    Equation e(node, firstNonZeroComp, [=](double time) -> double { return rhs(time) / direction[firstNonZeroComp]; });

    // add terms for all non zero direction components
    for (int iComponent = firstNonZeroComp + 1; iComponent < direction.rows(); ++iComponent)
        if (std::abs(direction[iComponent]) > 0)
            e.AddTerm(Term(node, iComponent, -direction[iComponent] / direction[firstNonZeroComp]));


    return e;
}

Equation Direction(const NodeSimple& node, Eigen::VectorXd direction, double value)
{
    return Direction(node, direction, RhsConstant(value));
}

std::vector<Equation> Direction(const Groups::Group<NodeSimple>& nodes, Eigen::VectorXd direction, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto& node : nodes)
        eqs.push_back(Direction(node, direction, rhs));
    return eqs;
}

std::vector<Equation> Direction(const Groups::Group<NodeSimple>& nodes, Eigen::VectorXd direction, double value)
{
    return Direction(nodes, direction, RhsConstant(value));
}

Equation Value(const NodeSimple& node, double value)
{
    return Value(node, [=](double) { return value; });
}

Equation Value(const NodeSimple& node, RhsFunction rhs)
{
    assert(node.GetNumValues() == 1 && "This function is ment to be used with single value nodes only");
    return Component(node, {eDirection::X}, rhs)[0];
}

std::vector<Equation> Value(const Groups::Group<NodeSimple>& nodes, double value)
{
    std::vector<Equation> eqs;
    for (auto& node : nodes)
        eqs.push_back(Value(node, value));
    return eqs;
}

std::vector<Equation> Value(const Groups::Group<NodeSimple>& nodes, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto& node : nodes)
        eqs.push_back(Value(node, rhs));
    return eqs;
}

} /* ConstraintPde */
} /* NuTo */
