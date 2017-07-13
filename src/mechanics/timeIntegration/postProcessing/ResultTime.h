#pragma once

#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

namespace NuTo
{
class ResultTime : public ResultBase
{
public:
    ResultTime(const std::string& rIdent);

    void CalculateAndAddValues(const StructureBase& rStructure, int timeStep,
                               const StructureOutputBlockVector& residual, double currentTime) override;

    //! @brief number of data points per time step (e.g. number of displacement components of a node)
    int GetNumData(const StructureBase& rStructure) const override
    {
        return 1;
    }

    std::unique_ptr<ResultBase> Clone() const override
    {
        return std::make_unique<ResultTime>(*this);
    }
};
} // namespace NuTo
