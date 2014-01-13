// $Id: Plane2D6N.cpp 276 2010-06-30 13:04:32Z arnold2 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <array>
#include <assert.h>
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/Plane2D6N.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::Plane2D6N::Plane2D6N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=6)
        throw MechanicsException("[NuTo::Plane2D6N::Plane2D6N] Exactly six nodes are required for this type of element.");
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
void NuTo::Plane2D6N::CalculateShapeFunctionsGeometry(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
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

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D6N::CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	CalculateShapeFunctionsGeometry(rNaturalCoordinates,rShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D6N::CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
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

    rDerivativeShapeFunctions[10] = -4.*s;
    rDerivativeShapeFunctions[11] = -8.*s-4.*r+4.;

}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D6N::CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	CalculateDerivativeShapeFunctionsGeometryNatural(rNaturalCoordinates,rDerivativeShapeFunctions);
}

//! @brief calculate the natural coordinates in 2D of all nodes
void NuTo::Plane2D6N::CalculateNaturalNodeCoordinates(std::vector< std::array<double,2> >& rNaturalNodeCoordinates)
{
	rNaturalNodeCoordinates.resize(6);
	//node 0
	(rNaturalNodeCoordinates[0])[0] = 0.0;
	rNaturalNodeCoordinates[0][1] = 0.0;

	//node 1
	rNaturalNodeCoordinates[1][0] = 1.0;
	rNaturalNodeCoordinates[1][1] = 0.0;

	//node 2
	rNaturalNodeCoordinates[2][0] = 0.0;
	rNaturalNodeCoordinates[2][1] = 1.0;

	//node 3
	rNaturalNodeCoordinates[3][0] = 0.5;
	rNaturalNodeCoordinates[3][1] = 0.0;

	//node 4
	rNaturalNodeCoordinates[4][0] = 0.5;
	rNaturalNodeCoordinates[4][1] = 0.5;

	//node 5
	rNaturalNodeCoordinates[5][0] = 0.0;
	rNaturalNodeCoordinates[5][1] = 0.5;
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Plane2D6N::CalculateShapeFunctionsSurface(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==3);
	rShapeFunctions[0] = 0.5*(1.-rNaturalCoordinates)-0.5*(1.-rNaturalCoordinates*rNaturalCoordinates);
	rShapeFunctions[1] = 1.-rNaturalCoordinates*rNaturalCoordinates;
	rShapeFunctions[2] = 0.5*(1.+rNaturalCoordinates)-0.5*(1.-rNaturalCoordinates*rNaturalCoordinates);
}

//! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D6N::CalculateDerivativeShapeFunctionsLocalSurface(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==3);
	rDerivativeShapeFunctions[0] = -0.5 + rNaturalCoordinates;
	rDerivativeShapeFunctions[1] = -2.0 * rNaturalCoordinates;
	rDerivativeShapeFunctions[2] =  0.5 + rNaturalCoordinates;
}

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Plane2D6N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{
	rSurfaceNodes.resize(3);
	switch(rSurface)
    {
    case 0:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[3];
    	rSurfaceNodes[2] = mNodes[1];
    	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[1];
    	rSurfaceNodes[1] = mNodes[4];
    	rSurfaceNodes[2] = mNodes[2];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[2];
    	rSurfaceNodes[1] = mNodes[5];
    	rSurfaceNodes[2] = mNodes[0];
    	break;
    default:
    	throw MechanicsException("[NuTo::Plane2D3N::GetSurfaceNodes] there are only 3 surfaces for a Plane2D3N.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Plane2D6N::GetNumSurfaces()const
{
    return 3;
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

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Plane2D6N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<6; count++)
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
template void NuTo::Plane2D6N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D6N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D6N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D6N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D6N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D6N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane2D6N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane2D6N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane2D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane2D6N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D6N)
#endif // ENABLE_SERIALIZATION
