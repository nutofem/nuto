// $Id: Plane2D15N.cpp 276 2010-06-30 13:04:32Z arnold2 $

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
#include "nuto/mechanics/elements/Plane2D15N.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::Plane2D15N::Plane2D15N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=15)
        throw MechanicsException("[NuTo::Plane2D15N::Plane2D15N] Exactly 10 nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    mNodes[3] = rNodes[3];
    mNodes[4] = rNodes[4];
    mNodes[5] = rNodes[5];
    mNodes[6] = rNodes[6];
    mNodes[7] = rNodes[7];
    mNodes[8] = rNodes[8];
    mNodes[9] = rNodes[9];
    mNodes[10] = rNodes[10];
    mNodes[11] = rNodes[11];
    mNodes[12] = rNodes[12];
    mNodes[13] = rNodes[13];
    mNodes[14] = rNodes[14];
    this->CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D15N::CalculateShapeFunctionsGeometry(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
    ShapeFunctions2D::ShapeFunctions2D15N(rNaturalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D15N::CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	CalculateShapeFunctionsGeometry(rNaturalCoordinates,rShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D15N::CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions2D::DerivativeShapeFunctions2D15N(rNaturalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D15N::CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	CalculateDerivativeShapeFunctionsGeometryNatural(rNaturalCoordinates,rDerivativeShapeFunctions);
}

//! @brief calculate the natural coordinates in 2D of all nodes
void NuTo::Plane2D15N::CalculateNaturalNodeCoordinates(std::vector< std::array<double,2> >& rNaturalNodeCoordinates)
{
	rNaturalNodeCoordinates.resize(15);

	std::array<double,2>& tmparray(rNaturalNodeCoordinates[0]);
	tmparray[0] = 0.;

	//node 0
	rNaturalNodeCoordinates[0][0] = 0.0;
	rNaturalNodeCoordinates[0][1] = 0.0;

	//node 1
	rNaturalNodeCoordinates[1][0] = 1./4.;
	rNaturalNodeCoordinates[1][1] = 0.0;

	//node 2
	rNaturalNodeCoordinates[2][0] = 2./4.;
	rNaturalNodeCoordinates[2][1] = 0.0;

	//node 3
	rNaturalNodeCoordinates[3][0] = 3./4.;
	rNaturalNodeCoordinates[3][1] = 0.0;

	//node 4
	rNaturalNodeCoordinates[4][0] = 1.;
	rNaturalNodeCoordinates[4][1] = 0.;

	//node 5
	rNaturalNodeCoordinates[5][0] = 0.;
	rNaturalNodeCoordinates[5][1] = 1./4.;

	//node 6
	rNaturalNodeCoordinates[6][0] = 1./4.;
	rNaturalNodeCoordinates[6][1] = 1./4.;

	//node 7
	rNaturalNodeCoordinates[7][0] = 2./4.;
	rNaturalNodeCoordinates[7][1] = 1./4.;

	//node 8
	rNaturalNodeCoordinates[8][0] = 3./4.;
	rNaturalNodeCoordinates[8][1] = 1./4.;

	//node 9
	rNaturalNodeCoordinates[9][0] = 0.;
	rNaturalNodeCoordinates[9][1] = 2./4.;

	//node 10
	rNaturalNodeCoordinates[10][0] = 1./4.;
	rNaturalNodeCoordinates[10][1] = 2./4.;

	//node 11
	rNaturalNodeCoordinates[11][0] = 2./4.;
	rNaturalNodeCoordinates[11][1] = 2./4.;

	//node 12
	rNaturalNodeCoordinates[12][0] = 0.;
	rNaturalNodeCoordinates[12][1] = 3./4.;

	//node 13
	rNaturalNodeCoordinates[13][0] = 1./4.;
	rNaturalNodeCoordinates[13][1] = 3./4.;

	//node 14
	rNaturalNodeCoordinates[14][0] = 0.;
	rNaturalNodeCoordinates[14][1] = 1.;
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Plane2D15N::CalculateShapeFunctionsSurface(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)const
{
    ShapeFunctions1D::ShapeFunctions1D5N(rNaturalCoordinates, rShapeFunctions);
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D15N::CalculateDerivativeShapeFunctionsLocalSurface(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions1D::DerivativeShapeFunctions1D5N(rNaturalCoordinates, rDerivativeShapeFunctions);
}

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Plane2D15N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{
	rSurfaceNodes.resize(5);
	switch(rSurface)
    {
    case 0:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[1];
    	rSurfaceNodes[2] = mNodes[2];
    	rSurfaceNodes[3] = mNodes[3];
    	rSurfaceNodes[4] = mNodes[4];
    	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[4];
    	rSurfaceNodes[1] = mNodes[8];
    	rSurfaceNodes[2] = mNodes[11];
    	rSurfaceNodes[3] = mNodes[13];
    	rSurfaceNodes[4] = mNodes[14];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[14];
    	rSurfaceNodes[1] = mNodes[12];
    	rSurfaceNodes[2] = mNodes[9];
    	rSurfaceNodes[3] = mNodes[5];
    	rSurfaceNodes[4] = mNodes[0];
    	break;
    default:
    	throw MechanicsException("[NuTo::Plane2D3N::GetSurfaceNodes] there are only 3 surfaces for a Plane2D3N.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Plane2D15N::GetNumSurfaces()const
{
    return 3;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Plane2D15N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType2D3NGauss13Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Plane2D15N::ReorderNodes()
{
	throw MechanicsException("[NuTo::Plane2D15N::ReorderNodes] not implemented.");
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Plane2D15N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<15; count++)
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
template void NuTo::Plane2D15N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D15N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D15N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D15N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D15N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D15N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane2D15N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane2D15N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane2D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane2D15N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D15N)
#endif // ENABLE_SERIALIZATION
