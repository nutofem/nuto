// $Id: LoadLoadSurfaceBase2D.cpp 178 2009-12-11 20:53:12Z eckardt4 $
#include <set>
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadSurfaceBase2D.h"
#include "nuto/mechanics/elements/Plane2D.h"


//! @brief constructor
NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D(int rLoadCase, StructureBase* rStructure, int rElementGroupId, int rNodeGroupId) : LoadBase(rLoadCase)
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

    std::cout << "number of loaded elements " << nodeGroup->GetNumMembers() << std::endl;

	//loop over all elements
	std::vector<const NodeBase*> surfaceNodes;
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	//check if plane element
        	Plane2D* elementPtr = itElement->second->AsPlane2D();

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
        			mPlaneElements.push_back (std::pair<const Plane2D*, int>(elementPtr,countSurface));
            		double nodeCoordinates[2];
            		surfaceNodes[0]->GetCoordinates2D(nodeCoordinates);
            		surfaceNodes[1]->GetCoordinates2D(nodeCoordinates);
        		}
        	}
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D] Error calculating surfaces for surface loads in element "  + ss.str() + "(Maybe not a solid element?).");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
        }
    }

    //set standard integration types for triangles and quads, this can be modified according to the needs
    mIntegrationType2NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss1Ip);
    mIntegrationType3NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss2Ip);
    mIntegrationType4NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss2Ip);
    mIntegrationType5NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss3Ip);
}

//! @brief adds the load to global sub-vectors
//! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
//! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
void NuTo::LoadSurfaceBase2D::AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rActiceDofsLoadVector, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const
{
    if (rLoadCase!=mLoadCase)
    	return;
	for (unsigned int countPlaneElement=0; countPlaneElement<mPlaneElements.size(); countPlaneElement++)
	{
		const Plane2D* planeElementPtr = mPlaneElements[countPlaneElement].first;

		std::vector<const NodeBase*> surfaceNodes;

		planeElementPtr->GetSurfaceNodes(mPlaneElements[countPlaneElement].second, surfaceNodes);
		int numSurfaceNodes(surfaceNodes.size());
		IntegrationTypeBase* integrationType(0);
		switch (numSurfaceNodes)
		{
		case 2:
			integrationType = mIntegrationType2NPtr;
			break;
		case 3:
			integrationType = mIntegrationType3NPtr;
			break;
		case 4:
			integrationType = mIntegrationType4NPtr;
			break;
		case 5:
			integrationType = mIntegrationType5NPtr;
			break;
		default:
			throw MechanicsException("[NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D] integration types only for 2, 3, 4 and 5 nodes (on the surface) implemented.");
		}

		//loop over surface integration points
        std::vector<double> shapeFunctions(numSurfaceNodes),derivativeShapeFunctionsLocal(numSurfaceNodes);
        for (int countIp=0; countIp<integrationType->GetNumIntegrationPoints(); countIp++)
		{
			//get local ip coordinates
        	double naturalIpCoordinates;
        	integrationType->GetLocalIntegrationPointCoordinates1D(countIp,naturalIpCoordinates);
			//std::cout << "naturalIpCoordinates " << naturalIpCoordinates << std::endl;

			//calculate shape functions
        	planeElementPtr->CalculateShapeFunctionsSurface(naturalIpCoordinates, shapeFunctions);
			//std::cout << "  shape functions " << std::endl;
			//for (unsigned int count=0; count<shapeFunctions.size(); count++)
			//	std::cout << shapeFunctions[count] << " ";
			//std::cout << std::endl;

			//calculate derivatives of shape functions
        	planeElementPtr->CalculateDerivativeShapeFunctionsLocalSurface(naturalIpCoordinates, derivativeShapeFunctionsLocal);

			//calculate global coordinates and basis vectors in the plane including the jacobian
        	FullVector<double,2> globalIpCoordinates,basisVector1, basisVector2;
        	globalIpCoordinates.setZero();
        	basisVector1.setZero();
        	basisVector2.setZero();
        	for (int countNode=0; countNode<numSurfaceNodes; countNode++)
        	{
        		double nodeCoordinates[2];
        		surfaceNodes[countNode]->GetCoordinates2D(nodeCoordinates);
    			//std::cout << "nodeCoordinates " << nodeCoordinates[0] << " " << nodeCoordinates[1] << std::endl;
        		globalIpCoordinates[0]+=shapeFunctions[countNode]*nodeCoordinates[0];
        		globalIpCoordinates[1]+=shapeFunctions[countNode]*nodeCoordinates[1];

				basisVector1(0) +=derivativeShapeFunctionsLocal[countNode]*nodeCoordinates[0];
				basisVector1(1) +=derivativeShapeFunctionsLocal[countNode]*nodeCoordinates[1];
        	}
			//std::cout << "globalIpCoordinates " << globalIpCoordinates[0] << " " << globalIpCoordinates[1]<< std::endl;

			//calculate basis vectors (in the plane) and determinant of Jacobian using the cross product of the basis vectors
        	basisVector2(0) = basisVector1(1);
        	basisVector2(1) = -basisVector1(0);
        	double detJ = basisVector1.norm();
        	basisVector2*=1./detJ;

        	//calculate weighting factor
        	double factor(detJ*(integrationType->GetIntegrationPointWeight(countIp)));

			//calculate surface load
        	FullVector<double,2> loadVector;
			CalculateSurfaceLoad(globalIpCoordinates, basisVector2, loadVector);
			loadVector*=factor;

			//std::cout << "load vector \n" << loadVector << std::endl;
			//std::cout << "  detJ " << detJ << " weight " << integrationType->GetIntegrationPointWeight(countIp) << std::endl;
			//std::cout << "  shape functions " << std::endl;
			//for (unsigned int count=0; count<shapeFunctions.size(); count++)
			//	std::cout << shapeFunctions[count] << " ";
			//std::cout << std::endl;

			//add load vector to global vector
			for (int countNode=0; countNode<numSurfaceNodes; countNode++)
			{
				assert(surfaceNodes[countNode]->GetNumDisplacements() == 2);
				for (int countDispDof=0; countDispDof<2; countDispDof++)
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
