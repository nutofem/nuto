#pragma once

#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

namespace NuTo
{
class StructureBase;
class NodeBase;

class ResultNodeDof : public ResultBase
{
public:
    ResultNodeDof(const std::string& rIdent, int rNodeId);

    //! @brief calculate the relevant nodal dofs and add to the internal routine
    void CalculateAndAddValues(const StructureBase& rStructure, int timeStep,
                               const StructureOutputBlockVector& residual, double currentTime) override;

    //! @brief calculate the relevant nodal dofs
    virtual Eigen::VectorXd CalculateValues(const StructureBase& rStructure) const = 0;

    void Info() const override;

protected:
    int mNodeId;
};
} // namespace NuTo
