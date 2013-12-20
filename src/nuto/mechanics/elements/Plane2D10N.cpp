// $Id: Plane2D10N.cpp 276 2010-06-30 13:04:32Z arnold2 $

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
#include "nuto/mechanics/elements/Plane2D10N.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::Plane2D10N::Plane2D10N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Plane2D::Plane2D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=10)
        throw MechanicsException("[NuTo::Plane2D10N::Plane2D10N] Exactly 10 nodes are required for this type of element.");
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
    this->CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D10N::CalculateShapeFunctionsGeometry(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==10);
    double r(rNaturalCoordinates[0]);
    double s(rNaturalCoordinates[1]);

    rShapeFunctions[0] =  + 1.0 -6.76014311053 * r -6.76014311053 * s + 11.5202862211 * r*r + 25.5608586632 * r*s + 11.5202862211 * s*s -5.76014311053 * r*r*r -19.8007155526 * r*r*s -19.8007155526 * r*s*s -5.76014311053 * s*s*s;
    rShapeFunctions[1] =  + 8.09016994375 * r -18.5103667207 * r*r -21.940482998 * r*s + 10.420196777 * r*r*r + 24.2705098312 * r*r*s + 13.8503130543 * r*s*s;
    rShapeFunctions[2] =  -2.33002683322 * r + 12.7502236102 * r*r -1.10008944408 * r*s -10.420196777 * r*r*r -6.99008049967 * r*r*s + 3.4301162773 * r*s*s;
    rShapeFunctions[3] =  + 1.0 * r -5.76014311053 * r*r + 2.52028622105 * r*s + 5.76014311053 * r*r*r -2.52028622105 * r*r*s -2.52028622105 * r*s*s;
    rShapeFunctions[4] =  + 8.09016994375 * s -21.940482998 * r*s -18.5103667207 * s*s + 13.8503130543 * r*r*s + 24.2705098312 * r*s*s + 10.420196777 * s*s*s;
    rShapeFunctions[5] =  + 27.0 * r*s -27.0 * r*r*s -27.0 * r*s*s;
    rShapeFunctions[6] =  -5.76014311053 * r*s + 13.8503130543 * r*r*s + 3.4301162773 * r*s*s;
    rShapeFunctions[7] =  -2.33002683322 * s -1.10008944408 * r*s + 12.7502236102 * s*s + 3.4301162773 * r*r*s -6.99008049967 * r*s*s -10.420196777 * s*s*s;
    rShapeFunctions[8] =  -5.76014311053 * r*s + 3.4301162773 * r*r*s + 13.8503130543 * r*s*s;
    rShapeFunctions[9] =  + 1.0 * s + 2.52028622105 * r*s -5.76014311053 * s*s -2.52028622105 * r*r*s -2.52028622105 * r*s*s + 5.76014311053 * s*s*s;
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane2D10N::CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	CalculateShapeFunctionsGeometry(rNaturalCoordinates,rShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D10N::CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==20);
    double r(rNaturalCoordinates[0]);
    double s(rNaturalCoordinates[1]);

    rDerivativeShapeFunctions[0] = -6.76014311053 + 23.0405724421*r + 25.5608586632*s-17.2804293316*r*r-39.6014311053*r*s-19.8007155526*s*s;
    rDerivativeShapeFunctions[1] = -6.76014311053 + 25.5608586632*r + 23.0405724421*s-19.8007155526*r*r-39.6014311053*r*s-17.2804293316*s*s;
    rDerivativeShapeFunctions[2] =  + 8.09016994375-37.0207334414*r-21.940482998*s + 31.2605903309*r*r + 48.5410196625*r*s + 13.8503130543*s*s;
    rDerivativeShapeFunctions[3] = -21.940482998*r + 24.2705098312*r*r + 27.7006261086*r*s;
    rDerivativeShapeFunctions[4] = -2.33002683322 + 25.5004472204*r-1.10008944408*s-31.2605903309*r*r-13.9801609993*r*s + 3.4301162773*s*s;
    rDerivativeShapeFunctions[5] = -1.10008944408*r-6.99008049967*r*r + 6.8602325546*r*s;
    rDerivativeShapeFunctions[6] =  + 1.0-11.5202862211*r + 2.52028622105*s + 17.2804293316*r*r-5.0405724421*r*s-2.52028622105*s*s;
    rDerivativeShapeFunctions[7] =  + 2.52028622105*r-2.52028622105*r*r-5.0405724421*r*s;
    rDerivativeShapeFunctions[8] = -21.940482998*s + 27.7006261086*r*s + 24.2705098312*s*s;
    rDerivativeShapeFunctions[9] =  + 8.09016994375-21.940482998*r-37.0207334414*s + 13.8503130543*r*r + 48.5410196625*r*s + 31.2605903309*s*s;
    rDerivativeShapeFunctions[10] =  + 27.0*s-54.0*r*s-27.0*s*s;
    rDerivativeShapeFunctions[11] =  + 27.0*r-27.0*r*r-54.0*r*s;
    rDerivativeShapeFunctions[12] = -5.76014311053*s + 27.7006261086*r*s + 3.4301162773*s*s;
    rDerivativeShapeFunctions[13] = -5.76014311053*r + 13.8503130543*r*r + 6.8602325546*r*s;
    rDerivativeShapeFunctions[14] = -1.10008944408*s + 6.8602325546*r*s-6.99008049967*s*s;
    rDerivativeShapeFunctions[15] = -2.33002683322-1.10008944408*r + 25.5004472204*s + 3.4301162773*r*r-13.9801609993*r*s-31.2605903309*s*s;
    rDerivativeShapeFunctions[16] = -5.76014311053*s + 6.8602325546*r*s + 13.8503130543*s*s;
    rDerivativeShapeFunctions[17] = -5.76014311053*r + 3.4301162773*r*r + 27.7006261086*r*s;
    rDerivativeShapeFunctions[18] =  + 2.52028622105*s-5.0405724421*r*s-2.52028622105*s*s;
    rDerivativeShapeFunctions[19] =  + 1.0 + 2.52028622105*r-11.5202862211*s-2.52028622105*r*r-5.0405724421*r*s + 17.2804293316*s*s;
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D10N::CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	CalculateDerivativeShapeFunctionsGeometryNatural(rNaturalCoordinates,rDerivativeShapeFunctions);
}

//! @brief calculate the natural coordinates in 2D of all nodes
void NuTo::Plane2D10N::CalculateNaturalNodeCoordinates(std::vector< std::array<double,2> >& rNaturalNodeCoordinates)
{
	rNaturalNodeCoordinates.resize(10);

	std::array<double,2>& tmparray(rNaturalNodeCoordinates[0]);
	tmparray[0] = 0.;

	//node 0
	rNaturalNodeCoordinates[0][0] = 0.0;
	rNaturalNodeCoordinates[0][1] = 0.0;

	//node 1
	rNaturalNodeCoordinates[1][0] = sqrt(1./20.);
	rNaturalNodeCoordinates[1][1] = 0.0;

	//node 2
	rNaturalNodeCoordinates[2][0] = 1.-sqrt(1./20.);
	rNaturalNodeCoordinates[2][1] = 0.0;

	//node 3
	rNaturalNodeCoordinates[3][0] = 1.;
	rNaturalNodeCoordinates[3][1] = 0.0;

	//node 4
	rNaturalNodeCoordinates[4][0] = 0.;
	rNaturalNodeCoordinates[4][1] = sqrt(1./20.);

	//node 5
	rNaturalNodeCoordinates[5][0] = 1./3.;
	rNaturalNodeCoordinates[5][1] = 1./3.;

	//node 6
	rNaturalNodeCoordinates[6][0] = 1.-sqrt(1./20.);
	rNaturalNodeCoordinates[6][1] = sqrt(1./20.);

	//node 7
	rNaturalNodeCoordinates[7][0] = 0;
	rNaturalNodeCoordinates[7][1] = 1.-sqrt(1./20.);

	//node 8
	rNaturalNodeCoordinates[8][0] = sqrt(1./20.);
	rNaturalNodeCoordinates[8][1] = 1.-sqrt(1./20.);

	//node 9
	rNaturalNodeCoordinates[9][0] = 0.;
	rNaturalNodeCoordinates[9][1] = 1.;
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Plane2D10N::CalculateShapeFunctionsSurface(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==4);
	double r(rNaturalCoordinates);
    double r2(rNaturalCoordinates*rNaturalCoordinates);
    double r3(r2*rNaturalCoordinates);
    rShapeFunctions[0] =  -0.220017888816 + 0.220017888816*r + 0.720017888816*r2 -0.720017888816*r3;
    rShapeFunctions[1] =  + 0.720017888816 -1.30252459712*r -0.720017888816*r2 + 1.30252459712*r3;
    rShapeFunctions[2] =  + 0.720017888816 + 1.30252459712*r -0.720017888816*r2 -1.30252459712*r3;
    rShapeFunctions[3] =  -0.220017888816 -0.220017888816*r + 0.720017888816*r2 + 0.720017888816*r3;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Plane2D10N::CalculateDerivativeShapeFunctionsLocalSurface(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==4);
	double r(rNaturalCoordinates);
    double r2(rNaturalCoordinates*rNaturalCoordinates);
    rDerivativeShapeFunctions[0] =  + 0.220017888816 + 1.44003577763*r-2.16005366645*r2;
    rDerivativeShapeFunctions[1] = -1.30252459712-1.44003577763*r + 3.90757379136*r2;
    rDerivativeShapeFunctions[2] =  + 1.30252459712-1.44003577763*r-3.90757379136*r2;
    rDerivativeShapeFunctions[3] = -0.220017888816 + 1.44003577763*r + 2.16005366645*r2;                                                                        }

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Plane2D10N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{
	rSurfaceNodes.resize(4);
	switch(rSurface)
    {
    case 0:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[1];
    	rSurfaceNodes[2] = mNodes[2];
    	rSurfaceNodes[3] = mNodes[3];
    	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[3];
    	rSurfaceNodes[1] = mNodes[6];
    	rSurfaceNodes[2] = mNodes[8];
    	rSurfaceNodes[3] = mNodes[9];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[9];
    	rSurfaceNodes[1] = mNodes[7];
    	rSurfaceNodes[2] = mNodes[4];
    	rSurfaceNodes[3] = mNodes[0];
    	break;
    default:
    	throw MechanicsException("[NuTo::Plane2D3N::GetSurfaceNodes] there are only 3 surfaces for a Plane2D3N.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Plane2D10N::GetNumSurfaces()const
{
    return 3;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Plane2D10N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType2D3NGauss13Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Plane2D10N::ReorderNodes()
{
	throw MechanicsException("[NuTo::Plane2D10N::ReorderNodes] not implemented.");
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Plane2D10N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<10; count++)
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
template void NuTo::Plane2D10N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D10N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D10N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane2D10N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D10N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane2D10N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane2D10N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane2D10N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane2D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane2D10N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D10N)
#endif // ENABLE_SERIALIZATION
