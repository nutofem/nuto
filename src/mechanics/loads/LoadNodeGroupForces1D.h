#pragma once

#include "mechanics/loads/LoadNodeGroup.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

namespace NuTo
{
template <class T>
class Group;

//! @brief Class for all forces applied to a group of nodes in 1D
class LoadNodeGroupForces1D : public LoadNodeGroup
{

public:
    //! @brief Constructor
    //! @param direction Direction of the force
    //! @param value Value of the force
    LoadNodeGroupForces1D(const Group<NodeBase>* rGroup, double direction, double value);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

protected:
    double mValue; //!< prescribed load of the node
    double mDirection; //!< direction of the force
};
} // namespace NuTo
