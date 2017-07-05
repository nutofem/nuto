#pragma once

#include "mechanics/loads/LoadNode.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

namespace NuTo
{

class LoadNodeForces2D : public LoadNode
{

public:
    //! @brief Constructor
    //! @param Irection Direction of the force
    //! @param value Value of the force
    LoadNodeForces2D(const NodeBase* node, const Eigen::MatrixXd& direction, double value);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

protected:
    double mValue;  //!< prescribed absolute value of the force at the node
    double mDirection[2]; //!< direction of the applied force (normalized)

private:
    LoadNodeForces2D() = default;
};
}//namespace NuTo

