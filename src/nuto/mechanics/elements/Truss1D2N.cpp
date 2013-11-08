// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/Truss1D2N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>

NuTo::Truss1D2N::Truss1D2N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Truss1D::Truss1D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
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
void NuTo::Truss1D2N::CalculateShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==2);
    rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates);
    rShapeFunctions[1] = 0.5*(1.+rLocalCoordinates);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D2N::CalculateShapeFunctionsField(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==2);
    rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates);
    rShapeFunctions[1] = 0.5*(1.+rLocalCoordinates);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D2N::CalculateShapeFunctionsNonlocalTotalStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==2);
    rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates);
    rShapeFunctions[1] = 0.5*(1.+rLocalCoordinates);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D2N::CalculateDerivativeShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==2);
    rDerivativeShapeFunctions[0] = -0.5;
    rDerivativeShapeFunctions[1] = 0.5;
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D2N::CalculateDerivativeShapeFunctionsField(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==2);
    rDerivativeShapeFunctions[0] = -0.5;
    rDerivativeShapeFunctions[1] = 0.5;
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D2N::CalculateDerivativeShapeFunctionsNonlocalTotalStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==2);
    rDerivativeShapeFunctions[0] = -0.5;
    rDerivativeShapeFunctions[1] = 0.5;
}


//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Truss1D2N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType1D2NGauss1Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Truss1D2N::ReorderNodes()
{
    std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[0];
    this->mNodes[0] = this->mNodes[1];
    this->mNodes[1] = tmp;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Truss1D2N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<2; count++)
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
template void NuTo::Truss1D2N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D2N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D2N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D2N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D2N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D2N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Truss1D2N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Truss1D2N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Truss1D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Truss1D2N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Truss1D2N)
#endif // ENABLE_SERIALIZATION
