// $Id: Plane2D4NSpectralOrder2.cpp 647 2013-10-31 13:23:00Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <tuple>

#include <assert.h>

#include <boost/foreach.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/Plane2D4NSpectralOrder2.h"
#include "nuto/mechanics/nodes/NodeBase.h"
NuTo::Plane2D4NSpectralOrder2::Plane2D4NSpectralOrder2(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=9)
        throw MechanicsException("[NuTo::Plane2D4NSpectralOrder2::Plane2D4NSpectralOrder2] Exactly 9 nodes are required for this type of element.");
    for (int count=0; count<9; count++)
    	mNodes[count] = rNodes[count];
    this->CheckElement(); //Geometry

    //check the position of the nodes on the edges
    std::cout << "check node position on edges" << std::endl;

}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D4NSpectralOrder2::CalculateShapeFunctionsGeometry(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==4);
    rShapeFunctions[0] = 0.25*(1.-rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
    rShapeFunctions[1] = 0.25*(1.+rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
    rShapeFunctions[2] = 0.25*(1.+rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
    rShapeFunctions[3] = 0.25*(1.-rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D4NSpectralOrder2::CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==9);
	double x1 = 0.5*rNaturalCoordinates[0]*(1.+rNaturalCoordinates[0]);
	double x2 = 1.-(rNaturalCoordinates[0]*rNaturalCoordinates[0]);
	double x3 = 0.5*rNaturalCoordinates[0]*(-1.+rNaturalCoordinates[0]);

	double y1 = 0.5*rNaturalCoordinates[1]*(1.+rNaturalCoordinates[1]);
	double y2 = 1.-(rNaturalCoordinates[1]*rNaturalCoordinates[1]);
	double y3 = 0.5*rNaturalCoordinates[1]*(-1.+rNaturalCoordinates[1]);

	rShapeFunctions[0] = x1*y1;
	rShapeFunctions[1] = x2*y1;
	rShapeFunctions[2] = x3*y1;
	rShapeFunctions[3] = x1*y2;
	rShapeFunctions[4] = x2*y2;
	rShapeFunctions[5] = x3*y2;
	rShapeFunctions[6] = x1*y3;
	rShapeFunctions[7] = x2*y3;
	rShapeFunctions[8] = x3*y3;
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D4NSpectralOrder2::CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
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

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D4NSpectralOrder2::CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==18);
	double x1 = 0.5*rNaturalCoordinates[0]*(1.+rNaturalCoordinates[0]);
	double x2 = 1.-(rNaturalCoordinates[0]*rNaturalCoordinates[0]);
	double x3 = 0.5*rNaturalCoordinates[0]*(-1.+rNaturalCoordinates[0]);

	double y1 = 0.5*rNaturalCoordinates[1]*(1.+rNaturalCoordinates[1]);
	double y2 = 1.-(rNaturalCoordinates[1]*rNaturalCoordinates[1]);
	double y3 = 0.5*rNaturalCoordinates[1]*(-1.+rNaturalCoordinates[1]);

	double dx1 = 0.5+rNaturalCoordinates[0];
	double dx2 = -2.*rNaturalCoordinates[0];
	double dx3 = -0.5+rNaturalCoordinates[0];

	double dy1 = 0.5+rNaturalCoordinates[1];
	double dy2 = -2.*rNaturalCoordinates[1];
	double dy3 = -0.5+rNaturalCoordinates[1];

    rDerivativeShapeFunctions[0] = dx1*y1;
    rDerivativeShapeFunctions[1] = x1*dy1;

    rDerivativeShapeFunctions[2] = dx2*y1;
    rDerivativeShapeFunctions[3] = x2*dy1;

    rDerivativeShapeFunctions[4] = dx3*y1;
    rDerivativeShapeFunctions[5] = x3*dy1;

    rDerivativeShapeFunctions[6] = dx1*y2;
    rDerivativeShapeFunctions[7] = x1*dy2;

    rDerivativeShapeFunctions[8] = dx2*y2;
    rDerivativeShapeFunctions[9] = x2*dy2;

    rDerivativeShapeFunctions[10] = dx3*y2;
    rDerivativeShapeFunctions[11] = x3*dy2;

    rDerivativeShapeFunctions[12] = dx1*y3;
    rDerivativeShapeFunctions[13] = x1*dy3;

    rDerivativeShapeFunctions[14] = dx2*y3;
    rDerivativeShapeFunctions[15] = x2*dy3;

    rDerivativeShapeFunctions[16] = dx3*y3;
    rDerivativeShapeFunctions[17] = x3*dy3;
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Plane2D4NSpectralOrder2::CalculateShapeFunctionsSurface(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==3);
    rShapeFunctions[0] = 0.5*rNaturalCoordinates*(1.+rNaturalCoordinates);
    rShapeFunctions[1] = 1.-(rNaturalCoordinates*rNaturalCoordinates);
    rShapeFunctions[2] = 0.5*rNaturalCoordinates*(-1.+rNaturalCoordinates);
}

//! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D4NSpectralOrder2::CalculateDerivativeShapeFunctionsLocalSurface(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==3);
    rDerivativeShapeFunctions[0] = 0.5+rNaturalCoordinates;
    rDerivativeShapeFunctions[1] = -2.*rNaturalCoordinates;
    rDerivativeShapeFunctions[2] = -0.5+rNaturalCoordinates;
}

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Plane2D4NSpectralOrder2::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{
	rSurfaceNodes.resize(3);
	switch(rSurface)
    {
    case 0:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[1];
    	rSurfaceNodes[2] = mNodes[2];
    	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[2];
    	rSurfaceNodes[1] = mNodes[5];
    	rSurfaceNodes[2] = mNodes[8];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[8];
    	rSurfaceNodes[1] = mNodes[7];
    	rSurfaceNodes[2] = mNodes[6];
    	break;
    case 3:
    	rSurfaceNodes[0] = mNodes[6];
    	rSurfaceNodes[1] = mNodes[3];
    	rSurfaceNodes[2] = mNodes[0];
    	break;
    default:
    	throw MechanicsException("[NuTo::Plane2D4NSpectralOrder2::GetSurfaceNodes] there are only 4 surfaces for a Plane2D4NSpectralOrder2.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Plane2D4NSpectralOrder2::GetNumSurfaces()const
{
    return 4;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Plane2D4NSpectralOrder2::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType2D4NLobatto9Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Plane2D4NSpectralOrder2::ReorderNodes()
{
    throw MechanicsException("[NuTo::Plane2D4NSpectralOrder2::ReorderNodes] not implemented.");
	//std::cout << "reorder element nodes" << std::endl;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Plane2D4NSpectralOrder2::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<9; count++)
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
template void NuTo::Plane2D4NSpectralOrder2::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D4NSpectralOrder2::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D4NSpectralOrder2::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D4NSpectralOrder2::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D4NSpectralOrder2::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D4NSpectralOrder2::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane2D4NSpectralOrder2::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane2D4NSpectralOrder2" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane2D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane2D4NSpectralOrder2" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D4NSpectralOrder2)
#endif // ENABLE_SERIALIZATION
