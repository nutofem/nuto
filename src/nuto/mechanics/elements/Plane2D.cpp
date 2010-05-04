// $Id: $

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/elements/Plane2D.h"
#include "nuto/mechanics/nodes/NodeCoordinates.h"
#include "nuto/mechanics/nodes/NodeDisplacements.h"
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
        const NodeCoordinates<2> *nodePtr(dynamic_cast<const NodeCoordinates<2> *>(GetNode(theNode)));
        assert(nodePtr!=0);
        nodePtr->GetCoordinates(&(rLocalCoordinates[2*theNode]));
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
        const NodeDisplacements<2> *nodePtr(dynamic_cast<const NodeDisplacements<2> *>(GetNode(theNode)));
        assert(nodePtr!=0);
        nodePtr->GetDisplacements(&(rLocalDisplacements[2*theNode]));
    }
}

// build global row dofs
void NuTo::Plane2D::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const
{
	rGlobalRowDofs.resize(2 * this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeDisplacementsBase *nodePtr = dynamic_cast< const NodeDisplacementsBase* >(this->GetNode(nodeCount));
        assert(nodePtr != NULL);
        rGlobalRowDofs[2 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
        rGlobalRowDofs[2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
    }
}

// build global column dof
void NuTo::Plane2D::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const
{
    this->CalculateGlobalRowDofs(rGlobalColumnDofs);
}

