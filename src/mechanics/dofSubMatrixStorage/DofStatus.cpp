/*
 * DofStatus.cpp
 *
 *  Created on: 4 Apr 2016
 *      Author: ttitsche
 */

#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/nodes/NodeEnum.h"
#include <boost/algorithm/string.hpp>


NuTo::DofStatus::DofStatus()
    : mHasInteractingConstraints(false)
{
}

namespace NuTo
{
std::ostream& operator<<(std::ostream& out, const DofStatus& dofStatus)
{
    out << "[---DofStatus:---]\n";
    out << "Existing DOF types:\n";
    for (auto dof : dofStatus.mDofTypes)
    {
        out << Node::DofToString(dof) << "; ";
    }
    out << "Active DOF types:\n";
    for (auto dof : dofStatus.mActiveDofTypes)
    {
        out << Node::DofToString(dof) << "; ";
    }
    for (auto dof : dofStatus.mNumActiveDofs)
    {
        out << "Number of active Dofs of type " << Node::DofToString(dof.first) << ": " << dof.second << "\n";
    }
    for (auto dof : dofStatus.mNumDependentDofs)
    {
        out << "Number of dependent Dofs of type " << Node::DofToString(dof.first) << ": " << dof.second << "\n";
    }
    return out;
}
} // namespace NuTo
