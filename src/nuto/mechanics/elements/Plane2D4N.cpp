// $Id$

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
#include "nuto/mechanics/elements/Plane2D4N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
NuTo::Plane2D4N::Plane2D4N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=4)
        throw MechanicsException("[NuTo::Plane2D4N::Plane2D4N] Exactly four nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    mNodes[3] = rNodes[3];
    this->CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D4N::CalculateShapeFunctionsGeometry(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
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
void NuTo::Plane2D4N::CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	CalculateShapeFunctionsGeometry(rNaturalCoordinates,rShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D4N::CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
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
void NuTo::Plane2D4N::CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	CalculateDerivativeShapeFunctionsGeometryNatural(rNaturalCoordinates,rDerivativeShapeFunctions);
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Plane2D4N::CalculateShapeFunctionsSurface(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==2);
    rShapeFunctions[0] = 0.5*(1.-rNaturalCoordinates);
    rShapeFunctions[1] = 0.5*(1.+rNaturalCoordinates);
}

//! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D4N::CalculateDerivativeShapeFunctionsLocalSurface(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==2);
    rDerivativeShapeFunctions[0] = -0.5;
    rDerivativeShapeFunctions[1] = 0.5;
}

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Plane2D4N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{
	rSurfaceNodes.resize(2);
	switch(rSurface)
    {
    case 0:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[1];
    	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[1];
    	rSurfaceNodes[1] = mNodes[2];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[2];
    	rSurfaceNodes[1] = mNodes[3];
    	break;
    case 3:
    	rSurfaceNodes[0] = mNodes[3];
    	rSurfaceNodes[1] = mNodes[0];
    	break;
    default:
    	throw MechanicsException("[NuTo::Plane2D4N::GetSurfaceNodes] there are only 4 surfaces for a Plane2D4N.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Plane2D4N::GetNumSurfaces()const
{
    return 4;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Plane2D4N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType2D4NGauss4Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Plane2D4N::ReorderNodes()
{
	//std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[3];
    this->mNodes[3] = tmp;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Plane2D4N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
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

//! @brief returns the natural coordinates of an given point
//! @param rGlobCoords (input) ... pointer to the array of coordinates
//! @param rLocCoords (output) ... coordinates to be returned
//! @return True if coordinates are within the element, False otherwise
bool NuTo::Plane2D4N::GetLocalPointCoordinates(const double rGlobCoords[2],  double rLocCoords[2])const
{
	//! local variables
	double inc[2] = {0.0,0.0};
	double iterPoint[3];
	double diffCoords[2] = {0.0,0.0};
	const double accuracy=1e-14;
	const unsigned int maxIter=100;
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumNodesGeometry());

	//! initialize local coordinates (starting value: center of element)
	rLocCoords[0] = 0.0;
	rLocCoords[1] = 0.0;

	//! get nodal coordinates
    std::vector<double> localNodeCoord(2*GetNumNodesGeometry());
    CalculateLocalCoordinates(localNodeCoord);

	//! check if node is inside this element
    //! if point outside the surrounding polygon: return false
    if( !CheckPointInside(rGlobCoords) )
    	return false;


	//! iterative solution of equations
	for(unsigned int i=0; ; ++i)
	{
		//! get inverse of jacobian
		double invJacobian[4], detJacobian;

	    InterpolateCoordinatesFrom2D( rLocCoords , iterPoint);
		diffCoords[0]=iterPoint[0] - rGlobCoords[0];
		diffCoords[1]=iterPoint[1] - rGlobCoords[1];

		// abort if accuracy is reached
		if( (diffCoords[0]*diffCoords[0]+diffCoords[1]*diffCoords[1]) < accuracy )
			return true;

		//! get inverse jacobian matrix
        CalculateDerivativeShapeFunctionsGeometryNatural(iterPoint, derivativeShapeFunctionsNatural);
        try
        {
        	CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, invJacobian, detJacobian);
        }
        catch ( MechanicsException  &e )
        {
        	throw MechanicsException("[NuTo::Plane2D4N::GetLocalPointCoordinates] Can't find the natural coordinates of the given point:  Determinant of the Jacobian is zero!");
        }

		//! jacInverse*(X_i-X_P)
		inc[0]=invJacobian[0]*diffCoords[0]+invJacobian[2]*diffCoords[1];
		inc[1]=invJacobian[1]*diffCoords[0]+invJacobian[3]*diffCoords[1];

		rLocCoords[0]-=inc[0];
		rLocCoords[1]-=inc[1];

		// error if no convergence within maximum number of iterations
		if(i==maxIter)
	        throw MechanicsException("[NuTo::Plane2D4N::GetLocalPointCoordinates] Can't find the natural coordinates of the given point: No convergence!.");
	}

	//! output
	//! if function had passed, point is outside the element -> FALSE will returned
	return false;
}

//! @brief checks if a node is inside of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! @param rGlobCoords (input) ... pointer to the array of coordinates
//! @return True if coordinates are within the element, False otherwise
bool NuTo::Plane2D4N::CheckPointInside(const double rGlobCoords[2])const
{
	const std::tuple<double,double> point(std::make_tuple(rGlobCoords[0],rGlobCoords[1]));
	std::vector<std::tuple<double,double> > polygon;
	BOOST_FOREACH(NuTo::NodeBase* thisNode,this->mNodes)
		polygon.push_back(std::make_tuple( thisNode->GetCoordinate(0),thisNode->GetCoordinate(1)));

	return CheckPointInsidePolygon(&point,&polygon);
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Plane2D4N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D4N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D4N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D4N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D4N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D4N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane2D4N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane2D4N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane2D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane2D4N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D4N)
#endif // ENABLE_SERIALIZATION
