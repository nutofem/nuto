/*
 * ResultDispNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "nuto/mechanics/timeIntegration/ResultNodeDisp.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::ResultNodeDisp::ResultNodeDisp(const std::string& rIdent, const NodeBase* rNodePtr) : ResultNodeDof(rIdent, rNodePtr)
{
}

//! @brief calculate the relevant nodal dofs
void NuTo::ResultNodeDisp::CalculateValues(NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)
{
    switch(mNodePtr->GetNumDisplacements())
    {
    case 1:
    {
    	mNodePtr->GetDisplacements1D(rValues.data());
    }
    break;
    case 2:
    {
    	mNodePtr->GetDisplacements2D(rValues.data());
    }
    break;
    case 3:
    {
    	mNodePtr->GetDisplacements3D(rValues.data());
    }
    break;
    default:
        break;
    }
}

//! @brief number of dofs (e.g. number of displacement components of a node
int NuTo::ResultNodeDisp::GetNumDofs()const
{
	return mNodePtr->GetNumDisplacements();
}


