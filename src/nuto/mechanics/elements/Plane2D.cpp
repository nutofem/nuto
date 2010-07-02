// $Id: $

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/elements/Plane2D.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>


NuTo::Plane2D::Plane2D(NuTo::StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        Plane(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{}

//! @brief calculates the local coordinates of the nodes
//! @param localCoordinates vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Plane2D::CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const
{
	assert((int)rLocalCoordinates.size()==2*GetNumNodes());
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        GetNode(theNode)->GetCoordinates2D(&(rLocalCoordinates[2*theNode]));
    }
}

//! @brief calculates the local displacements of the nodes
//! @param localDisplacements vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Plane2D::CalculateLocalDisplacements(std::vector<double>& rLocalDisplacements)const
{
	assert((int)rLocalDisplacements.size()==2*GetNumNodes());
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        GetNode(theNode)->GetDisplacements2D(&(rLocalDisplacements[2*theNode]));
    }
}

// build global row dofs
void NuTo::Plane2D::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const
{
	rGlobalRowDofs.resize(2 * this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeBase *nodePtr = this->GetNode(nodeCount);
        rGlobalRowDofs[2 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
        rGlobalRowDofs[2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
    }
}

// build global column dof
void NuTo::Plane2D::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const
{
    int NumNonlocalElements(GetNumNonlocalElements());
	if (GetNumNonlocalElements()==0)
	    this->CalculateGlobalRowDofs(rGlobalColumnDofs);
    else
    {
	    const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
        int NumCols(0);
        for (int theNonlocalElement=0; theNonlocalElement<NumNonlocalElements; theNonlocalElement++)
        {
        	NumCols += nonlocalElements[theNonlocalElement]->AsPlane()->GetNumLocalDofs();
        }
        rGlobalColumnDofs.resize(NumCols);
        int shift(0);
        for (int theNonlocalElement=0; theNonlocalElement<NumNonlocalElements; theNonlocalElement++)
        {
            const ElementBase* nonlocalElement(nonlocalElements[theNonlocalElement]);
        	for (int nodeCount = 0; nodeCount < nonlocalElement->GetNumNodes(); nodeCount++)
            {
                const NodeBase *nodePtr = nonlocalElement->GetNode(nodeCount);
                rGlobalColumnDofs[shift + 2 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
                rGlobalColumnDofs[shift + 2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
            }
            shift+=2*nonlocalElement->GetNumNodes();
        }
        assert(shift==NumCols);
    }
}

//! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
//! @return Area
double NuTo::Plane2D::CalculateArea()const
{
    double coordinates1[2];
    double coordinates2[2];
    GetNode(GetNumNodes()-1)->GetCoordinates2D(coordinates2);
	double area(0);
	for (int theNode = 0; theNode<GetNumNodes(); theNode++)
	{
		coordinates1[0] = coordinates2[0];
		coordinates1[1] = coordinates2[1];
		GetNode(theNode)->GetCoordinates2D(coordinates2);
		area += coordinates1[0]*coordinates2[1] - coordinates1[1]*coordinates2[0];
	}
    return 0.5*area;
}


