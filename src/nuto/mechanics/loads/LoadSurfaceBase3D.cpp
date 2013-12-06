// $Id: LoadLoadSurfaceBase3D.cpp 178 2009-12-11 20:53:12Z eckardt4 $

#include <set>
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadSurfaceBase3D.h"
#include "nuto/mechanics/elements/Solid.h"


//! @brief constructor
NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D(int rLoadCase, StructureBase* rStructure, int rElementGroupId, int rNodeGroupId) : LoadBase(rLoadCase)
{
    //get element group
	const Group<ElementBase> *elementGroup = rStructure->GroupGetGroupPtr(rElementGroupId)->AsGroupElement();

    //get node group
	const Group<NodeBase> *nodeGroup = rStructure->GroupGetGroupPtr(rNodeGroupId)->AsGroupNode();

	//since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
	std::set<const NodeBase*> nodePtrSet;
    for (Group<NodeBase>::const_iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
    	nodePtrSet.insert(itNode->second);
    }

	//loop over all elements
	std::vector<const NodeBase*> surfaceNodes;
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	//check if solid element
        	Solid* elementPtr = itElement->second->AsSolid();

        	//loop over all surfaces
        	for (int countSurface=0; countSurface<elementPtr->GetNumSurfaces(); countSurface++)
        	{
        		bool addSurface(true);
        		elementPtr->GetSurfaceNodes(countSurface, surfaceNodes);

        		//check, if all surface nodes are in the node group
        		for (unsigned int countNode=0; countNode<surfaceNodes.size(); countNode++)
        		{
        			if (nodePtrSet.find(surfaceNodes[countNode])==nodePtrSet.end())
        			{
        				//this surface has at least on node that is not in the list, continue
        				addSurface=false;
        			}
        		}

        		if (addSurface)
        		{
        			mVolumeElements.push_back (std::pair<const Solid*, int>(elementPtr,countSurface));
        		}
        	}
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D] Error calculating surfaces for surface loads in element "  + ss.str() + "(Maybe not a solid element?).");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
        }

    }

    //set standard integration types for triangles and quads, this can be modified according to the needs
    mIntegrationType3NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType2D3NGauss1Ip);
    mIntegrationType4NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType2D4NGauss4Ip);
    mIntegrationType6NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType2D3NGauss3Ip);
}

//! @brief adds the load to global sub-vectors
//! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
//! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
void NuTo::LoadSurfaceBase3D::AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rActiceDofsLoadVector, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
	for (unsigned int countVolumeElement=0; countVolumeElement<mVolumeElements.size(); countVolumeElement++)
	{
		const Solid* solidElementPtr = mVolumeElements[countVolumeElement].first;

		std::vector<const NodeBase*> surfaceNodes;

		solidElementPtr->GetSurfaceNodes(mVolumeElements[countVolumeElement].second, surfaceNodes);
		int numSurfaceNodes(surfaceNodes.size());
		IntegrationTypeBase* integrationType(0);
		switch (numSurfaceNodes)
		{
		case 3:
			integrationType = mIntegrationType3NPtr;
			break;
		case 4:
			integrationType = mIntegrationType4NPtr;
			break;
		case 6:
			integrationType = mIntegrationType6NPtr;
			break;
		default:
			throw MechanicsException("[NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D] integration types only for 3 4 and 6 nodes implemented.");
		}

		//loop over surface integration points
        std::vector<double> shapeFunctions(numSurfaceNodes),derivativeShapeFunctionsLocal(2*numSurfaceNodes);
        for (int countIp=0; countIp<integrationType->GetNumIntegrationPoints(); countIp++)
		{
			//get local ip coordinates
        	double localIpCoordinates[2];
        	integrationType->GetLocalIntegrationPointCoordinates2D(countIp,localIpCoordinates);

			//calculate shape functions
			solidElementPtr->CalculateShapeFunctionsSurface(localIpCoordinates, shapeFunctions);

			//calculate derivatives of shape functions
			solidElementPtr->CalculateDerivativeShapeFunctionsLocalSurface(localIpCoordinates, derivativeShapeFunctionsLocal);

			//calculate global coordinates and basis vectors in the plane including the jacobian
        	FullVector<double,3> globalIpCoordinates, basisVector1, basisVector2, basisVector3;
        	basisVector1.setZero();
        	basisVector2.setZero();
        	for (int countNode=0; countNode<numSurfaceNodes; countNode++)
        	{
        		double nodeCoordinates[3];
        		surfaceNodes[countNode]->GetCoordinates3D(nodeCoordinates);
				globalIpCoordinates[0]+=shapeFunctions[countNode]*nodeCoordinates[0];
				globalIpCoordinates[1]+=shapeFunctions[countNode]*nodeCoordinates[1];
				globalIpCoordinates[2]+=shapeFunctions[countNode]*nodeCoordinates[2];

				basisVector1(0) +=derivativeShapeFunctionsLocal[2*countNode]*nodeCoordinates[0];
				basisVector1(1) +=derivativeShapeFunctionsLocal[2*countNode]*nodeCoordinates[1];
				basisVector1(2) +=derivativeShapeFunctionsLocal[2*countNode]*nodeCoordinates[2];

				basisVector2(0) +=derivativeShapeFunctionsLocal[2*countNode+1]*nodeCoordinates[0];
				basisVector2(1) +=derivativeShapeFunctionsLocal[2*countNode+1]*nodeCoordinates[1];
				basisVector2(2) +=derivativeShapeFunctionsLocal[2*countNode+1]*nodeCoordinates[2];
        	}

			//calculate basis vectors (in the plane) and determinant of Jacobian using the cross product of the basis vectors
        	basisVector3(0) = basisVector1(1)*basisVector2(2)-basisVector1(2)*basisVector2(1);
        	basisVector3(1) = basisVector1(2)*basisVector2(0)-basisVector1(0)*basisVector2(2);
        	basisVector3(2) = basisVector1(0)*basisVector2(1)-basisVector1(1)*basisVector2(0);
        	double detJ = basisVector3.norm();
        	basisVector3*=1./detJ;

        	//calculate weighting factor
        	double factor(detJ*(integrationType->GetIntegrationPointWeight(countIp)));

			//calculate surface load
        	FullVector<double,3> loadVector;
			CalculateSurfaceLoad(globalIpCoordinates, basisVector3, loadVector);
			loadVector*=factor;

			//std::cout << "load vector \n" << loadVector << std::endl;
			//std::cout << "  detJ " << detJ << " weight " << integrationType->GetIntegrationPointWeight(countIp) << std::endl;
			//std::cout << "  shape functions " << std::endl;
			//for (int count=0; count<shapeFunctions.size(); count++)
			//	std::cout << shapeFunctions[count] << " ";
			//std::cout << std::endl;

			//add load vector to global vector
			for (int countNode=0; countNode<numSurfaceNodes; countNode++)
			{
				assert(surfaceNodes[countNode]->GetNumDisplacements() == 3);
				for (int countDispDof=0; countDispDof<3; countDispDof++)
				{
					int theDof = surfaceNodes[countNode]->GetDofDisplacement(countDispDof);
					if (theDof<rActiceDofsLoadVector.GetNumRows())
					{
						//std::cout << "add to dof " << theDof << " " << shapeFunctions[countNode]*loadVector(countDispDof) << std::endl;
						rActiceDofsLoadVector(theDof)+=shapeFunctions[countNode]*loadVector(countDispDof);
					}
					else
					{
						rDependentDofsLoadVector(theDof-rActiceDofsLoadVector.GetNumRows())+=shapeFunctions[countNode]*loadVector(countDispDof);
					}
				}
			}
		}
	}
}
