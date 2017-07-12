#pragma once

#include "mechanics/timeIntegration/postProcessing/ResultNodeDof.h"

namespace NuTo
{
class ResultNodeDisp : public ResultNodeDof
{
public:
    ResultNodeDisp(const std::string& rIdent, int rNodeId);

    //! @brief calculate the relevant nodal dofs
    Eigen::VectorXd CalculateValues(const StructureBase& rStructure) const override;

    //! @brief number of data points per time step (e.g. number of displacement components of a node
    int GetNumData(const StructureBase& rStructure) const override;

    NuTo::eTimeIntegrationResultType GetResultType() const override;

    std::string GetTypeId() const
    {
        return std::string("ResultNodeDisp");
    }

    void Info() const override
    {
    }

    std::unique_ptr<ResultBase> Clone() const override
    {
        return std::make_unique<ResultNodeDisp>(*this);
    }
};
}
