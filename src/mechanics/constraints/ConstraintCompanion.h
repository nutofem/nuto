#pragma once

#include "mechanics/constraints/Constraints.h"
#include "mechanics/groups/GroupBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/DirectionEnum.h"

namespace NuTo
{
namespace Constraint
{

RhsFunction RhsConstant(double value)
{
    return [=](double) { return value; };
}

RhsFunction RhsRamp(double timeEnd, double valueEnd)
{
    return [=](double time) { return valueEnd * time / timeEnd; };
}

std::vector<Equation> Fix(const NodeBase& node, std::vector<eDirection> directions)
{
    std::vector<Equation> eqs;
    for (auto direction : directions)
    {
        int component = ToComponentIndex(direction);
        if (node.GetNumDofs() < component)
            // TODO This check is not meaningful at the moment, since this method returns the
            // total number of dofs, say 4 (3 disp, 1 temp). This will not find the error
            // if you try to constrain the Z component of the temperature. Which would be wrong.
            throw MechanicsException(__PRETTY_FUNCTION__, "Dimension mismatch");

        eqs.push_back(Equation(RhsConstant(0.), {Term(node, component, 1)}));
    }
    return eqs;
}

std::vector<Equation> Fix(const Group<NodeBase>& nodes, std::vector<eDirection> directions)
{
    std::vector<Equation> eqs;
    for (auto& nodePair : nodes)
    {
        auto tmpEqs = Fix(*nodePair.second, directions);
        eqs.insert(eqs.end(), tmpEqs.begin(), tmpEqs.end());
    }
    return eqs;
}

Equation Fix(const NodeBase& node)
{
    return Equation(RhsConstant(0.), {Term(node, 0, 1)});
}

std::vector<Equation> Fix(const Group<NodeBase>& nodes)
{
    std::vector<Equation> eqs;
    for (auto& nodePair : nodes)
        eqs.push_back(Fix(*nodePair.second));
    return eqs;
}

Equation Direction(const NodeBase& node, Eigen::VectorXd direction, RhsFunction rhs)
{
    if (node.GetNumDofs() != direction.rows())
        // TODO This check is not meaningful at the moment, since this method returns the
        // total number of dofs, say 4 (3 disp, 1 temp). This will not find the error
        // if you try to constrain the Z component of the temperature. Which would be wrong.
        throw MechanicsException(__PRETTY_FUNCTION__, "Dimension mismatch");

    direction.normalize();
    Equation e(rhs);
    for (int iComponent = 0; iComponent < direction.rows(); ++iComponent)
        e.AddTerm(Term(node, iComponent, direction[iComponent]));

    return e;
}

std::vector<Equation> Direction(const Group<NodeBase>& nodes, Eigen::VectorXd direction, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto& nodePair : nodes)
        eqs.push_back(Direction(*nodePair.second, direction, rhs));
    return eqs;
}

Equation Value(const NodeBase& node, RhsFunction rhs)
{
    return Direction(node, Eigen::VectorXd::Ones(1), rhs);
}

std::vector<Equation> Value(const Group<NodeBase>& nodes, RhsFunction rhs)
{
    return Direction(nodes, Eigen::VectorXd::Ones(1), rhs);
}

} /* Constraint */
} /* NuTo */
