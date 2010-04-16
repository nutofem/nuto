// $Id: $

#include "nuto/mechanics/elements/Truss1D2N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include <assert.h>

NuTo::Truss1D2N::Truss1D2N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes, ElementDataBase::eElementDataType rElementDataType) :
        NuTo::Truss1D::Truss1D(rStructure, rElementDataType, GetStandardIntegrationType())
{
    if (rNodes.size()!=2)
        throw MechanicsException("[NuTo::Truss2N::Truss2N] Exactly two nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    this->CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D2N::CalculateShapeFunctions(const double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==2);
    rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates);
    rShapeFunctions[1] = 0.5*(1.+rLocalCoordinates);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D2N::CalculateDerivativeShapeFunctions(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==2);
    rDerivativeShapeFunctions[0] = -0.5;
    rDerivativeShapeFunctions[1] = 0.5;
}


//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationTypeBase::eIntegrationType NuTo::Truss1D2N::GetStandardIntegrationType()
{
    return NuTo::IntegrationTypeBase::IntegrationType1D2NGauss1Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Truss1D2N::ReorderNodes()
{
    std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[0];
    this->mNodes[0] = this->mNodes[1];
    this->mNodes[1] = tmp;
}
