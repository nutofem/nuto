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
    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot);

    //! @brief calculate the relevant nodal dofs
    virtual Eigen::VectorXd CalculateValues(const StructureBase& rStructure) const = 0;

    ResultNodeDof* AsResultNodeDof() override
    {
        return this;
    }

    void Info() const override;

protected:
    int mNodeId;
};
} // namespace NuTo
