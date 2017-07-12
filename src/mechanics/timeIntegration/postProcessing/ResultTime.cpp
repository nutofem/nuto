#include "mechanics/timeIntegration/postProcessing/ResultTime.h"

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

using namespace NuTo;

ResultTime::ResultTime(const std::string& rIdent)
    : ResultBase(rIdent)
{
}

void ResultTime::Info() const
{
}

void ResultTime::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot, double rTime)
{
    assert(rTimeStepPlot >= 0);
    if (rTimeStepPlot >= mData.rows())
    {
        this->Resize(rStructure, 2 * (rTimeStepPlot + 1), false);
    }
    mData(rTimeStepPlot, 0) = rTime;
}

eTimeIntegrationResultType ResultTime::GetResultType() const
{
    return eTimeIntegrationResultType::TIME;
}
