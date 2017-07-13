#pragma once

#include "mechanics/constraints/Constraints.h"
#include "mechanics/DirectionEnum.h"

namespace NuTo
{
template <typename T>
class Group;

namespace Constraint
{

//! @brief linear function
//! @param timeEnd ... end time
//! @param valueEnd ... end value
//! @return ... linear function crossing (0,0) and (timeEnd, valueEnd)
inline RhsFunction RhsRamp(double timeEnd, double valueEnd)
{
    return [=](double time) { return valueEnd * time / timeEnd; };
}


//! @brief constant function
//! @param constantValue ... constant value
//! @return ... constant function f(t) = constant value
inline RhsFunction RhsConstant(double constantValue)
{
    return [=](double time) { return constantValue; };
}


//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param node ... node reference
//! @param directions ... vector of directions (X, Y, Z)
//! @param value ... constant value
//! @return ... vector of constraint equations
std::vector<Equation> Component(const NodeBase& node, std::vector<eDirection> directions, double value = 0.0);

//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param nodes ... group of nodes
//! @param directions ... vector of directions (X, Y, Z)
//! @param value ... constant value
//! @return ... vector of constraint equations
std::vector<Equation> Component(const Group<NodeBase>& nodes, std::vector<eDirection> directions, double value = 0.0);

//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param node ... node reference
//! @param directions ... vector of directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Component(const NodeBase& node, std::vector<eDirection> directions, RhsFunction rhs);

//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param nodes ... group of nodes
//! @param directions ... vector of directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Component(const Group<NodeBase>& nodes, std::vector<eDirection> directions, RhsFunction rhs);

//! @param node ... node reference
//! @param direction ... directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... constraint equation
Equation Direction(const NodeBase& node, Eigen::VectorXd direction, RhsFunction rhs);

//! @param node ... node reference
//! @param direction ... directions (X, Y, Z)
//! @param value ... constant value
//! @return ... constraint equation
Equation Direction(const NodeBase& node, Eigen::VectorXd direction, double value = 0.0);

//! @param nodes ... group of nodes
//! @param direction ... directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Direction(const Group<NodeBase>& nodes, Eigen::VectorXd direction, RhsFunction rhs);

//! @param nodes ... group of nodes
//! @param direction ... directions (X, Y, Z)
//! @param value ... constant value
//! @return vector of constraint equations
std::vector<Equation> Direction(const Group<NodeBase>& nodes, Eigen::VectorXd direction, double value = 0.0);

//! @param node ... node reference
//! @param value ... constant value
//! @return ... constraint equation
Equation Value(const NodeBase& node, double value = 0.0);

//! @param node ... node reference
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... constraint equation
Equation Value(const NodeBase& node, RhsFunction rhs);

//! @param nodes ... group of nodes
//! @param value ... constant value
//! @return ... vector of constraint equations
std::vector<Equation> Value(const Group<NodeBase>& nodes, double value = 0.0);

//! @param nodes ... group of nodes
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Value(const Group<NodeBase>& nodes, RhsFunction rhs);

} /* Constraint */
} /* NuTo */
