#include "mechanics/timeIntegration/postProcessing/ResultTime.h"

#include "mechanics/nodes/NodeBase.h"

using namespace NuTo;

ResultTime::ResultTime(const std::string& rIdent)
    : ResultBase(rIdent)
{
}


void ResultTime::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot,
                                       const StructureOutputBlockVector& residual, double rTime)
{
    assert(rTimeStepPlot >= 0);
    if (rTimeStepPlot >= mData.rows())
    {
        Resize(rStructure, 2 * (rTimeStepPlot + 1), false);
    }
    mData(rTimeStepPlot, 0) = rTime;
}
