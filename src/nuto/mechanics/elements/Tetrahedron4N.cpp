#include "nuto/mechanics/elements/Tetrahedron4N.h"

NuTo::Tetrahedron4N::Tetrahedron4N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Solid::Solid(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
    if (rNodes.size()!=4)
        throw MechanicsException("[NuTo::Tetrahedron4N::Tetrahedron4N] Exactly four nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    mNodes[3] = rNodes[3];
}

NuTo::Tetrahedron4N::~Tetrahedron4N()
{}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Tetrahedron4N::CalculateShapeFunctions(const double rLocalCoordinates[3],
        std::vector <double> &rShapeFunctions) const
{
    assert(rShapeFunctions.size()==4);
    rShapeFunctions[0] = 1.0 - rLocalCoordinates[0]- rLocalCoordinates[1]
                         - rLocalCoordinates[2];
    rShapeFunctions[1] = rLocalCoordinates[0];
    rShapeFunctions[2] = rLocalCoordinates[1];
    rShapeFunctions[3] = rLocalCoordinates[2];
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Tetrahedron4N::CalculateDerivativeShapeFunctionsLocal(
    const double rLocalCoordinates[3],
    std::vector <double> &rDerivativeShapeFunctions) const
{
    assert(rDerivativeShapeFunctions.size()==12);

    //node1
    rDerivativeShapeFunctions[0] = -1.0;
    rDerivativeShapeFunctions[1] = -1.0;
    rDerivativeShapeFunctions[2] = -1.0;
    //node2
    rDerivativeShapeFunctions[3] = 1.0;
    rDerivativeShapeFunctions[4] = 0.0;
    rDerivativeShapeFunctions[5] = 0.0;
    //node3
    rDerivativeShapeFunctions[6] = 0.0;
    rDerivativeShapeFunctions[7] = 1.0;
    rDerivativeShapeFunctions[8] = 0.0;
    //node4
    rDerivativeShapeFunctions[9] = 0.0;
    rDerivativeShapeFunctions[10] = 0.0;
    rDerivativeShapeFunctions[11] = 1.0;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Tetrahedron4N::GetStandardIntegrationType()
{
    throw MechanicsException("Tetrahedron4N::getStandardIntegrationType: necessary integration type not implemented yet");
    return IntegrationType::IntegrationType3D8NGauss2x2x2Ip;
}

//! @brief reorder element nodes
void NuTo::Tetrahedron4N::ReorderNodes()
{
    // swap nodes 2 and 3
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[2];
    this->mNodes[2] = tmp;
}
