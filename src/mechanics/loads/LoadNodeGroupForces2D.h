#pragma once

#include "mechanics/loads/LoadNodeGroup.h"

namespace NuTo
{
class NodeBase;
template <class T>
class Group;
class StructureOutputBlockVector;

//! @brief Class for all forces applied to a group of nodes in 2D
class LoadNodeGroupForces2D : public LoadNodeGroup
{

public:
    //! @brief Constructor
    //! @param direction Direction of the force
    //! @param value Value of the force
    LoadNodeGroupForces2D(const Group<NodeBase>* group, const Eigen::VectorXd& direction, double value);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

protected:
    double mValue; //!< prescribed force of the node
    double mDirection[2]; //!< direction of the applied constraint
};
} // namespace NuTo
