#include "mechanics/timeIntegration/postProcessing/ResultNodeDisp.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

ResultNodeDisp::ResultNodeDisp(const std::string& rIdent, int rNodeId)
    : ResultNodeDof(rIdent, rNodeId)
{
}


Eigen::VectorXd ResultNodeDisp::CalculateValues(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->Get(Node::eDof::DISPLACEMENTS);
}


int ResultNodeDisp::GetNumData(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->GetNum(Node::eDof::DISPLACEMENTS);
}
