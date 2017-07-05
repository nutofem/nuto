/*
 * ResultDispNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */


#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/ResultNodeDisp.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

NuTo::ResultNodeDisp::ResultNodeDisp(const std::string& rIdent, int rNodeId)
    : ResultNodeDof(rIdent, rNodeId)
{
}

Eigen::VectorXd NuTo::ResultNodeDisp::CalculateValues(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->Get(Node::eDof::DISPLACEMENTS);
}

int NuTo::ResultNodeDisp::GetNumData(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->GetNum(Node::eDof::DISPLACEMENTS);
}

NuTo::eTimeIntegrationResultType NuTo::ResultNodeDisp::GetResultType() const
{
    return NuTo::eTimeIntegrationResultType::NODE_DISPLACEMENT;
}
