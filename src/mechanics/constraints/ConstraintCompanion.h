#pragma once

#include "mechanics/constraints/Constraints.h"
#include "mechanics/groups/Group.h"
#include "mechanics/DirectionEnum.h"

namespace NuTo
{
namespace Constraint
{

inline RhsFunction RhsRamp(double timeEnd, double valueEnd)
{
    return [=](double time) { return valueEnd * time / timeEnd; };
}

inline RhsFunction RhsConstant(double constantValue)
{
    return [=](double time) { return constantValue; }; 
}

std::vector<Equation> Component(const NodeBase& node, std::vector<eDirection> directions, double value = 0.0);

std::vector<Equation> Component(const Group<NodeBase>& nodes, std::vector<eDirection> directions, double value = 0.0);

std::vector<Equation> Component(const NodeBase& node, std::vector<eDirection> directions, RhsFunction rhs);

std::vector<Equation> Component(const Group<NodeBase>& nodes, std::vector<eDirection> directions,  RhsFunction rhs);

Equation Direction(const NodeBase& node, Eigen::VectorXd direction, RhsFunction rhs);

Equation Direction(const NodeBase& node, Eigen::VectorXd direction, double value = 0.0);

std::vector<Equation> Direction(const Group<NodeBase>& nodes, Eigen::VectorXd direction, RhsFunction rhs);

std::vector<Equation> Direction(const Group<NodeBase>& nodes, Eigen::VectorXd direction, double value = 0.0);

Equation Value(const NodeBase& node, double value = 0.0);

Equation Value(const NodeBase& node, RhsFunction rhs);

std::vector<Equation> Value(const Group<NodeBase>& nodes, double value = 0.0);

std::vector<Equation> Value(const Group<NodeBase>& nodes, RhsFunction rhs);

} /* Constraint */
} /* NuTo */
