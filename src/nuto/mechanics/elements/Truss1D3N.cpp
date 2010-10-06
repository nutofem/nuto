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
#include "nuto/mechanics/elements/Truss1D3N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include <assert.h>

NuTo::Truss1D3N::Truss1D3N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Truss1D::Truss1D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
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
void NuTo::Truss1D3N::CalculateShapeFunctions(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
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
void NuTo::Truss1D3N::CalculateDerivativeShapeFunctions(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==3);
	rDerivativeShapeFunctions[0] = -0.5 + rLocalCoordinates;
	rDerivativeShapeFunctions[1] = -2.0 * rLocalCoordinates;
	rDerivativeShapeFunctions[2] =  0.5 + rLocalCoordinates;
}


//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Truss1D3N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType1D2NGauss2Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Truss1D3N::ReorderNodes()
{
    std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[0];
    this->mNodes[0] = this->mNodes[2];
    this->mNodes[2] = tmp;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Truss1D3N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Truss1D3N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Truss1D3N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Truss1D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Truss1D3N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Truss1D3N)
#endif // ENABLE_SERIALIZATION
