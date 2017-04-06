#pragma once

#include "mechanics/constraints/Constraints.h"
#include "mechanics/groups/GroupBase.h"
#include "mechanics/groups/Group.h"

namespace NuTo
{
namespace Constraint
{


Equation FixVector(const NodeBase& node, Node::eDof dof, Eigen::VectorXd direction, RhsFunction rhs)
{
    int numComponents = node.GetNum(dof);
    if (numComponents != direction.rows())
        throw MechanicsException(__PRETTY_FUNCTION__, "Dimension mismatch");

    direction.normalize();
    Equation e(rhs);
    for (int iComponent = 0; iComponent < numComponents; ++iComponent)
        e.AddTerm(Term(node, iComponent, direction[iComponent]));

    return e;
}

Equation FixScalar(const NodeBase& node, Node::eDof dof, RhsFunction rhs)
{
    return FixVector(node, dof, Eigen::VectorXd::Ones(1), rhs);
}

//! @remark This should certainly be a Group<NodeBase>! However, this caused some wierd compilation problems, maybe
//! caused by foward declarations. So, this currently works and the groups need a cleanup anyways
std::vector<Equation> FixVector(const GroupBase& nodes, Node::eDof dof, Eigen::VectorXd direction, RhsFunction rhs)
{
    std::vector<Equation> equations;
    for (auto nodePair : *nodes.AsGroupNode())
        equations.push_back(FixVector(*nodePair.second, dof, direction, rhs));
    return equations;
}

std::vector<Equation> FixScalar(const GroupBase& nodes, Node::eDof dof, RhsFunction rhs)
{
    return FixVector(nodes, dof, Eigen::VectorXd::Ones(1), rhs);
}


} /* Constraint */
} /* NuTo */
