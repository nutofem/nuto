// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/Plane2D3N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include <assert.h>

NuTo::Plane2D3N::Plane2D3N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=3)
        throw MechanicsException("[NuTo::Plane2D3N::Plane2D3N] Exactly three nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    this->CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D3N::CalculateShapeFunctions(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==3);
    rShapeFunctions[0] = 1-rNaturalCoordinates[0]-rNaturalCoordinates[1];
    rShapeFunctions[1] = rNaturalCoordinates[0];
    rShapeFunctions[2] = rNaturalCoordinates[1];
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D3N::CalculateDerivativeShapeFunctionsNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==6);
    rDerivativeShapeFunctions[0] = -1.0;
    rDerivativeShapeFunctions[1] = -1.0;

    rDerivativeShapeFunctions[2] = 1.0;
    rDerivativeShapeFunctions[3] = 0.0;

    rDerivativeShapeFunctions[4] = 0.0;
    rDerivativeShapeFunctions[5] = 1.0;
}


//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Plane2D3N::GetStandardIntegrationType()
{
	return NuTo::IntegrationType::IntegrationType2D3NGauss1Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Plane2D3N::ReorderNodes()
{
	//std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[2];
    this->mNodes[2] = tmp;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Plane2D3N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<3; count++)
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
template void NuTo::Plane2D3N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D3N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D3N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D3N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D3N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D3N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane2D3N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane2D3N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane2D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane2D3N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D3N)
#endif // ENABLE_SERIALIZATION

