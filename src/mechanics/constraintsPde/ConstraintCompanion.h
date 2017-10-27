#pragma once

#include "base/Group.h"
#include "mechanics/constraintsPde/Constraints.h"
#include "mechanics/constraintsPde/Equation.h"
#include "mechanics/DirectionEnum.h"
#include "mechanics/nodes/NodeSimple.h"

namespace NuTo
{
template <typename T>
class Group;

namespace ConstraintPde
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
    return [=](double) { return constantValue; };
}


//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param node ... node reference
//! @param directions ... vector of directions (X, Y, Z)
//! @param value ... constant value
//! @return ... vector of constraint equations
std::vector<Equation> Component(const NodeSimple& node, std::vector<eDirection> directions, double value = 0.0);

//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param nodes ... group of nodes
//! @param directions ... vector of directions (X, Y, Z)
//! @param value ... constant value
//! @return ... vector of constraint equations
std::vector<Equation> Component(const Groups::Group<NodeSimple>& nodes, std::vector<eDirection> directions,
                                double value = 0.0);

//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param node ... node reference
//! @param directions ... vector of directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Component(const NodeSimple& node, std::vector<eDirection> directions, RhsFunction rhs);

//! @brief constraints components of vector valued nodes in X and/or Y and/or Z
//! @param nodes ... group of nodes
//! @param directions ... vector of directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Component(const Groups::Group<NodeSimple>& nodes, std::vector<eDirection> directions,
                                RhsFunction rhs);


//! @brief creates a constraint equation for constraints that are not axes aligned
//! @param node ... node reference
//! @param direction ... directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... constraint equation
Equation Direction(const NodeSimple& node, Eigen::VectorXd direction, RhsFunction rhs);

//! @brief creates a constraint equation for constraints that are not axes aligned
//! @param node ... node reference
//! @param direction ... directions (X, Y, Z)
//! @param value ... constant value
//! @return ... constraint equation
Equation Direction(const NodeSimple& node, Eigen::VectorXd direction, double value = 0.0);

//! @brief creates multiple constraint equations for constraints that are not axes aligned
//! @param nodes ... group of nodes
//! @param direction ... directions (X, Y, Z)
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Direction(const Groups::Group<NodeSimple>& nodes, Eigen::VectorXd direction, RhsFunction rhs);

//! @brief creates multiple constraint equations for constraints that are not axes aligned
//! @param nodes ... group of nodes
//! @param direction ... directions (X, Y, Z)
//! @param value ... constant value
//! @return vector of constraint equations
std::vector<Equation> Direction(const Groups::Group<NodeSimple>& nodes, Eigen::VectorXd direction, double value = 0.0);

//! @brief Constraint single value node
//! @param node ... node reference
//! @param value ... constant value
//! @return ... constraint equation
Equation Value(const NodeSimple& node, double value = 0.0);

//! @brief Constraint single value node
//! @param node ... node reference
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... constraint equation
Equation Value(const NodeSimple& node, RhsFunction rhs);

//! @brief Constraint group of single value nodes
//! @param nodes ... group of nodes
//! @param value ... constant value
//! @return ... vector of constraint equations
std::vector<Equation> Value(const Groups::Group<NodeSimple>& nodes, double value = 0.0);

//! @brief Constraint group of single value nodes
//! @param nodes ... group of nodes
//! @param rhs ... time dependent constraint function (double time) --> double
//! @return ... vector of constraint equations
std::vector<Equation> Value(const Groups::Group<NodeSimple>& nodes, RhsFunction rhs);

} /* Constraint */
} /* NuTo */
