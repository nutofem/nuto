#include "mechanics/timeIntegration/postProcessing/ResultNodeDof.h"

#include "mechanics/nodes/NodeBase.h"

using namespace NuTo;

ResultNodeDof::ResultNodeDof(const std::string& rIdent, int rNodeId)
    : ResultBase(rIdent)
{
    mNodeId = rNodeId;
}


void ResultNodeDof::Info() const
{
}


void ResultNodeDof::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot)
{
    assert(rTimeStepPlot >= 0);
    Eigen::VectorXd dofValues = this->CalculateValues(rStructure);
    if (rTimeStepPlot >= mData.rows())
    {
        this->Resize(rStructure, 2 * (rTimeStepPlot + 1), false);
    }
    if (dofValues.rows() != mData.cols())
        throw MechanicsException(__PRETTY_FUNCTION__, "the allocated number of rows is wrong.");
    mData.row(rTimeStepPlot) = dofValues.transpose();
}
