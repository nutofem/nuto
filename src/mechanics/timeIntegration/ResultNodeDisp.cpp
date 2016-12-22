/*
 * ResultDispNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "math/FullMatrix.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/ResultNodeDisp.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

NuTo::ResultNodeDisp::ResultNodeDisp(const std::string& rIdent, int rNodeId) : ResultNodeDof(rIdent, rNodeId)
{
}

//! @brief calculate the relevant nodal dofs
void NuTo::ResultNodeDisp::CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)const
{
    const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));

    rValues = node->Get(Node::eDof::DISPLACEMENTS).transpose();
}

//! @brief number of data points per time step (e.g. number of displacement components of a node
int NuTo::ResultNodeDisp::GetNumData(const StructureBase& rStructure)const
{
    const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));
    return node->GetNum(Node::eDof::DISPLACEMENTS);
}

NuTo::eTimeIntegrationResultType NuTo::ResultNodeDisp::GetResultType() const
{
    return NuTo::eTimeIntegrationResultType::NODE_DISPLACEMENT;
}


