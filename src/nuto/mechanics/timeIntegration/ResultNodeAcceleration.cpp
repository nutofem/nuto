/*
 * ResultAccelerationNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "nuto/mechanics/timeIntegration/ResultNodeAcceleration.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::ResultNodeAcceleration::ResultNodeAcceleration(const std::string& rIdent, int rNodeId) : ResultNodeDof(rIdent, rNodeId)
{
}

//! @brief calculate the relevant nodal dofs
void NuTo::ResultNodeAcceleration::CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)const
{
	const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));

	switch(node->GetNumDisplacements())
    {
    case 1:
    {
    	node->GetDisplacements1D(2,rValues.data());
    }
    break;
    case 2:
    {
    	node->GetDisplacements2D(2,rValues.data());
    }
    break;
    case 3:
    {
    	node->GetDisplacements3D(2,rValues.data());
    }
    break;
    default:
        break;
    }
}

//! @brief number of data points per time step (e.g. number of Accelerationlacement components of a node
int NuTo::ResultNodeAcceleration::GetNumData(const StructureBase& rStructure)const
{
	const NodeBase* node(rStructure.NodeGetNodePtr(mNodeId));
	return node->GetNumDisplacements();
}


