// $Id: $

#include "nuto/mechanics/elements/Plane2D4N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include <assert.h>

NuTo::Plane2D4N::Plane2D4N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes, ElementDataBase::eElementDataType rElementDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType())
{
	if (rNodes.size()!=4)
        throw MechanicsException("[NuTo::Plane2D4N::Plane2D4N] Exactly four nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    mNodes[3] = rNodes[3];
    this->CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D4N::CalculateShapeFunctions(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==4);
    rShapeFunctions[0] = 0.25*(1.-rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
    rShapeFunctions[1] = 0.25*(1.+rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
    rShapeFunctions[2] = 0.25*(1.+rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
    rShapeFunctions[3] = 0.25*(1.-rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D4N::CalculateDerivativeShapeFunctionsNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==8);
    rDerivativeShapeFunctions[0] = -0.25*(1.-rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[1] = -0.25*(1.-rNaturalCoordinates[0]);

    rDerivativeShapeFunctions[2] = +0.25*(1.-rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[3] = -0.25*(1.+rNaturalCoordinates[0]);

    rDerivativeShapeFunctions[4] = +0.25*(1.+rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[5] = +0.25*(1.+rNaturalCoordinates[0]);

    rDerivativeShapeFunctions[6] = -0.25*(1.+rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[7] = +0.25*(1.-rNaturalCoordinates[0]);
}


//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationTypeBase::eIntegrationType NuTo::Plane2D4N::GetStandardIntegrationType()
{
    return NuTo::IntegrationTypeBase::IntegrationType2D4NGauss4Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Plane2D4N::ReorderNodes()
{
	std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[3];
    this->mNodes[3] = tmp;
}
