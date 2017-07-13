/*
 * ResultAccelerationNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */


#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/ResultNodeAcceleration.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

NuTo::ResultNodeAcceleration::ResultNodeAcceleration(const std::string& rIdent, int rNodeId)
    : ResultNodeDof(rIdent, rNodeId)
{
}

Eigen::VectorXd NuTo::ResultNodeAcceleration::CalculateValues(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->Get(Node::eDof::DISPLACEMENTS, 2);
}

int NuTo::ResultNodeAcceleration::GetNumData(const StructureBase& rStructure) const
{
    return rStructure.NodeGetNodePtr(mNodeId)->GetNum(Node::eDof::DISPLACEMENTS);
}

NuTo::eTimeIntegrationResultType NuTo::ResultNodeAcceleration::GetResultType() const
{
    return NuTo::eTimeIntegrationResultType::NODE_ACCELERATION;
}
