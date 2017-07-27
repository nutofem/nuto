#pragma once

#include "mechanics/loads/LoadNode.h"

namespace NuTo
{
class StructureOutputBlockVector;

class LoadNodeForces3D : public LoadNode
{

public:
    //! @brief Constructor
    //! @param direction Direction of the force
    //! @param value Value of the force
    LoadNodeForces3D(const NodeBase* node, const Eigen::MatrixXd& direction, double value);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

protected:
    double mValue; //!< prescribed absolute value of the force at the node
    double mDirection[3]; //!< direction of the applied force (normalized)
};
} // namespace NuTo
