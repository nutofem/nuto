#pragma once

#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

namespace NuTo
{
class ResultTime : public ResultBase
{
public:
    ResultTime(const std::string& rIdent);

    std::string GetTypeId() const
    {
        return std::string("ResultTime");
    }

    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot, double rTime);

    NuTo::eTimeIntegrationResultType GetResultType() const override;

    //! @brief number of data points per time step (e.g. number of displacement components of a node)
    int GetNumData(const StructureBase& rStructure) const override
    {
        return 1;
    }

    ResultTime* AsResultTime() override
    {
        return this;
    }

    void Info() const override;

    std::unique_ptr<ResultBase> Clone() const override
    {
        return std::make_unique<ResultTime>(*this);
    } 
};
} // namespace NuTo
