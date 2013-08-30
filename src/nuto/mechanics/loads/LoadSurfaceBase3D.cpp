// $Id: LoadLoadSurfaceBase3D.cpp 178 2009-12-11 20:53:12Z eckardt4 $
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadSurfaceBase3D.h"
#include "nuto/mechanics/elements/Solid.h"


//! @brief constructor
NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D(int rLoadCase, const Group<ElementBase>* rElementGroup, const Group<NodeBase>* rNodeGroup) : LoadBase(rLoadCase)
{
	//determine the elements and corresponding surfaces by checking for all element surfaces if all nodes are in the node group
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

		//loop over surface integration points
        std::vector<double> shapeFunctions(numSurfaceNodes),derivativeShapeFunctionsLocal(2*numSurfaceNodes);
        for (int countIp=0; countIp<mIntegrationTypePtr->GetNumIntegrationPoints(); countIp++)
		{
			//get local ip coordinates
        	double localIpCoordinates[2];
        	mIntegrationTypePtr->GetLocalIntegrationPointCoordinates2D(countIp,localIpCoordinates);

			//calculate shape functions
			solidElementPtr->CalculateShapeFunctionsSurface(localIpCoordinates, shapeFunctions);

			//calculate derivatives of shape functions
			solidElementPtr->CalculateDerivativeShapeFunctionsLocalSurface(localIpCoordinates, derivativeShapeFunctionsLocal);

			//calculate global coordinates and basis vectors in the plane including the jacobian
        	FullVector<double,3> globalIpCoordinates, basisVector1, basisVector2, basisVector3;
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
        	basisVector3(0) = basisVector2(1)*basisVector3(2)-basisVector2(2)*basisVector3(1);
        	basisVector3(1) = basisVector2(2)*basisVector3(0)-basisVector2(0)*basisVector3(2);
        	basisVector3(2) = basisVector2(0)*basisVector3(1)-basisVector2(1)*basisVector3(0);
        	double detJ = basisVector3.norm();
        	basisVector3*=1./detJ;

        	//calculate weighting factor
        	double factor(detJ*(mIntegrationTypePtr->GetIntegrationPointWeight(countIp)));

			//calculate surface load
        	FullVector<double,3> loadVector;
			CalculateSurfaceLoad(globalIpCoordinates, basisVector3, loadVector);
			loadVector*=factor;

			//add load vector to global vector
			for (int countNode; countNode<numSurfaceNodes; countNode++)
			{
				assert(surfaceNodes[count]->GetNumDisplacements() == 3);
				for (int countDispDof=0; countDispDof<3; countDispDof++)
				{
					int theDof = surfaceNodes[countNode]->GetDofDisplacement(countDispDof);
					if (theDof<rActiceDofsLoadVector.GetNumRows())
					{
						rActiceDofsLoadVector(theDof)+=loadVector(countDispDof);
					}
					else
					{
						rDependentDofsLoadVector(theDof-rActiceDofsLoadVector.GetNumRows())+=loadVector(countDispDof);
					}
				}
			}
		}
	}
}
