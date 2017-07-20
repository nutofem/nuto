#pragma once

#include "mechanics/loads/LoadNode.h"

namespace NuTo
{
class StructureOutputBlockVector;

//! @brief Class for all forces applied to a single node in 1D
class LoadNodeForces1D : public LoadNode
{

public:
    //! @brief Constructor
    //! @param direction Direction of the force
    //! @param value Value of the force
    LoadNodeForces1D(const NodeBase* node, double direction, double value);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

protected:
    double mValue; //!< prescribed force of the node
    double mDirection; //!< direction of the force

private:
    LoadNodeForces1D() = default;
};
} // namespace NuTo
