#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeDof.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"

using namespace NuTo;

ResultGroupNodeDof::ResultGroupNodeDof(const std::string& rIdent, int rGroupNodeId)
    : ResultBase(rIdent)
{
    mGroupNodeId = rGroupNodeId;
}

void ResultGroupNodeDof::Info() const
{
}

void ResultGroupNodeDof::CalculateAndAddValues(const StructureBase& rStructure, int timeStep, const StructureOutputBlockVector& residual,
                           double currentTime)
{
    assert(timeStep >= 0);
    if (timeStep >= mData.rows())
    {
        Resize(rStructure, 2 * (timeStep + 1), false);
    }
    Eigen::VectorXd values = CalculateValues(rStructure, residual);

    if (values.rows() != mData.cols())
        throw MechanicsException(__PRETTY_FUNCTION__, "the allocated number of rows is wrong.");

    mData.row(timeStep) = values.transpose();
}

const Group<NodeBase>* ResultGroupNodeDof::GetGroupNodePtr(const StructureBase& rStructure) const
{
    return rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();
}
