// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/Tetrahedron4N.h"
#include "nuto/mechanics/nodes/NodeBase.h"

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
    CheckElement();
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

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Tetrahedron4N::CalculateShapeFunctionsSurface(const double rLocalCoordinates[2], std::vector<double>& rShapeFunctions)const
{

}

//! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Tetrahedron4N::CalculateDerivativeShapeFunctionsLocalSurface(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const
{

}

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Tetrahedron4N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Tetrahedron4N::GetNumSurfaces()const
{
    return 4;
}

//! @brief returns the order of the shape functions (required to determine the integration order for surface loads)
//! @return polynomial order
int NuTo::Tetrahedron4N::GetOrderOfShapeFunctions()const
{
    return 1;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Tetrahedron4N::GetStandardIntegrationType()
{
    return IntegrationType::IntegrationType3D4NGauss1Ip;
}

//! @brief reorder element nodes
void NuTo::Tetrahedron4N::ReorderNodes()
{
    // swap nodes 2 and 3
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[2];
    this->mNodes[2] = tmp;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Tetrahedron4N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<4; count++)
    {
        if (this->mNodes[count]==rOldPtr)
        {
            this->mNodes[count]=rNewPtr;
            break;
        }
    }
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Tetrahedron4N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron4N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron4N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron4N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron4N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron4N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Tetrahedron4N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Tetrahedron4N" << std::endl;
#endif
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Solid)
           & BOOST_SERIALIZATION_NVP(mNodes);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Tetrahedron4N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Tetrahedron4N)
#endif // ENABLE_SERIALIZATION
