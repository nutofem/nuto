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

NuTo::ResultNodeAcceleration::ResultNodeAcceleration(const std::string& rIdent, int rNodeId) : ResultNodeDof(rIdent, rNodeId)
{
}

//! @brief calculate the relevant nodal dofs
void NuTo::ResultNodeAcceleration::CalculateValues(const StructureBase& rStructure, Eigen::Matrix<double, 1, Eigen::Dynamic>& rValues)const
{
    const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));

    rValues = node->Get(Node::eDof::DISPLACEMENTS, 2).transpose();
}

//! @brief number of data points per time step (e.g. number of Accelerationlacement components of a node
int NuTo::ResultNodeAcceleration::GetNumData(const StructureBase& rStructure)const
{
    const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));
    return node->GetNum(Node::eDof::DISPLACEMENTS);
}

NuTo::eTimeIntegrationResultType NuTo::ResultNodeAcceleration::GetResultType() const
{
    return NuTo::eTimeIntegrationResultType::NODE_ACCELERATION;
}


