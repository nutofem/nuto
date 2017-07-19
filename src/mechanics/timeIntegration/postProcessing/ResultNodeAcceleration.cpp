#include "mechanics/timeIntegration/postProcessing/ResultNodeAcceleration.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

ResultNodeAcceleration::ResultNodeAcceleration(const std::string& rIdent, int rNodeId)
    : ResultNodeDof(rIdent, rNodeId)
{
}


Eigen::VectorXd ResultNodeAcceleration::CalculateValues(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->Get(Node::eDof::DISPLACEMENTS, 2);
}


int ResultNodeAcceleration::GetNumData(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->GetNum(Node::eDof::DISPLACEMENTS);
}
