// $Id: Plane2D3N.cpp 276 2010-06-30 13:04:32Z arnold2 $

#include "nuto/mechanics/elements/Plane2D6N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include <assert.h>

NuTo::Plane2D6N::Plane2D6N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=3)
        throw MechanicsException("[NuTo::Plane2D6N::Plane2D6N] Exactly three nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    mNodes[3] = rNodes[3];
    mNodes[4] = rNodes[4];
    mNodes[5] = rNodes[5];
    this->CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D6N::CalculateShapeFunctions(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==6);
    double r(rNaturalCoordinates[0]);
    double s(rNaturalCoordinates[1]);

    rShapeFunctions[0] =  2.*(r*r+s*s)+4.*r*s-3.*(r+s)+1.;
    rShapeFunctions[1] =  2.*r*r-r;
    rShapeFunctions[2] =  2.*s*s-s;
    rShapeFunctions[3] = -4.*r*(r+s-1.);
    rShapeFunctions[4] =  4.*r*s;
    rShapeFunctions[5] = -4.*s*(s+r-1.);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D6N::CalculateDerivativeShapeFunctionsNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==12);
    double r(rNaturalCoordinates[0]);
    double s(rNaturalCoordinates[1]);

    rDerivativeShapeFunctions[0] = 4.*(r+s)-3. ;
    rDerivativeShapeFunctions[1] = rDerivativeShapeFunctions[0];

    rDerivativeShapeFunctions[2] = 4.*r-1.;
    rDerivativeShapeFunctions[3] = 0.;

    rDerivativeShapeFunctions[4] = 0.;
    rDerivativeShapeFunctions[5] = 4.*s-1.;

    rDerivativeShapeFunctions[6] = -8.*r-4.*s+4.;
    rDerivativeShapeFunctions[7] = -4.*r;

    rDerivativeShapeFunctions[8] = 4.*s;
    rDerivativeShapeFunctions[9] = 4.*r;

    rDerivativeShapeFunctions[10] = -4.*s;;
    rDerivativeShapeFunctions[11] = -8.*s-4.*r+4.;

}


//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Plane2D6N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType2D3NGauss3Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Plane2D6N::ReorderNodes()
{
	//std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[2];
    this->mNodes[2] = tmp;
    tmp = this->mNodes[3];
    this->mNodes[3] = this->mNodes[5];
    this->mNodes[5] = tmp;

}
