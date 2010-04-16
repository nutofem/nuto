// $Id: $

#include "nuto/mechanics/elements/Truss1D3N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include <assert.h>

NuTo::Truss1D3N::Truss1D3N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes, ElementDataBase::eElementDataType rElementDataType) :
 NuTo::Truss1D::Truss1D(rStructure, rElementDataType, GetStandardIntegrationType())
{
	if (rNodes.size()!=3)
    	throw MechanicsException("[NuTo::Truss3N::Truss3N] Exactly three nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D3N::CalculateShapeFunctions(const double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==3);
	rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates)-0.5*(1.-rLocalCoordinates*rLocalCoordinates);
	rShapeFunctions[1] = 1.-rLocalCoordinates*rLocalCoordinates;
	rShapeFunctions[2] = 0.5*(1.+rLocalCoordinates)-0.5*(1.-rLocalCoordinates*rLocalCoordinates);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D3N::CalculateDerivativeShapeFunctions(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==3);
	rDerivativeShapeFunctions[0] = -0.5 + rLocalCoordinates;
	rDerivativeShapeFunctions[1] = -2.0 * rLocalCoordinates;
	rDerivativeShapeFunctions[2] =  0.5 + rLocalCoordinates;
}


//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationTypeBase::eIntegrationType NuTo::Truss1D3N::GetStandardIntegrationType()
{
    return NuTo::IntegrationTypeBase::IntegrationType1D2NGauss2Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Truss1D3N::ReorderNodes()
{
    std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[0];
    this->mNodes[0] = this->mNodes[2];
    this->mNodes[2] = tmp;
}
