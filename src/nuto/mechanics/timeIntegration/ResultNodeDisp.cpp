/*
 * ResultDispNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "nuto/mechanics/timeIntegration/ResultNodeDisp.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::ResultNodeDisp::ResultNodeDisp(const std::string& rIdent, int rNodeId) : ResultNodeDof(rIdent, rNodeId)
{
}

//! @brief calculate the relevant nodal dofs
void NuTo::ResultNodeDisp::CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)const
{
	const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));

	switch(node->GetNumDisplacements())
    {
    case 1:
    {
    	node->GetDisplacements1D(rValues.data());
    }
    break;
    case 2:
    {
    	node->GetDisplacements2D(rValues.data());
    }
    break;
    case 3:
    {
    	node->GetDisplacements3D(rValues.data());
    }
    break;
    default:
        break;
    }
}

//! @brief number of data points per time step (e.g. number of displacement components of a node
int NuTo::ResultNodeDisp::GetNumData(const StructureBase& rStructure)const
{
	const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));
	return node->GetNumDisplacements();
}


