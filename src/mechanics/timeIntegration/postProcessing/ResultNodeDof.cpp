#include "mechanics/timeIntegration/postProcessing/ResultNodeDof.h"

#include "mechanics/nodes/NodeBase.h"

using namespace NuTo;

ResultNodeDof::ResultNodeDof(const std::string& rIdent, int rNodeId)
    : ResultBase(rIdent)
{
    mNodeId = rNodeId;
}


void ResultNodeDof::CalculateAndAddValues(const StructureBase& structure, int timeStep,
                                          const StructureOutputBlockVector& residual, double currentTime)
{
    assert(timeStep >= 0);
    Eigen::VectorXd dofValues = CalculateValues(structure);
    if (timeStep >= mData.rows())
    {
        Resize(structure, 2 * (timeStep + 1), false);
    }
    if (dofValues.rows() != mData.cols())
        throw Exception(__PRETTY_FUNCTION__, "The allocated number of rows is wrong.");
    mData.row(timeStep) = dofValues.transpose();
}
