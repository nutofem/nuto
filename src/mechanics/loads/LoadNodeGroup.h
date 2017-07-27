#pragma once

#include "mechanics/loads/LoadBase.h"

namespace NuTo
{
template <class T>
class Group;
class NodeBase;

//! @brief Abstract class for all loads applied to a node group
class LoadNodeGroup : public LoadBase
{

public:
    //! @brief constructor
    LoadNodeGroup(const Group<NodeBase>* rGroup);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumConstraintEquations() const;


protected:
    const Group<NodeBase>* mGroup;
};
} // namespace NuTo
