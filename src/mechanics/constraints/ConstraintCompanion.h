#pragma once

#include "mechanics/constraints/Constraints.h"

namespace NuTo
{
namespace Constraint
{

Equation FixDof(const NodeBase& node, Node::eDof dof, Eigen::VectorXd direction, double rhs)
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



} /* Constraint */
} /* NuTo */
