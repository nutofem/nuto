// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/Brick8N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include <assert.h>

NuTo::Brick8N::Brick8N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Solid::Solid(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
    if (rNodes.size()!=8)
        throw MechanicsException("[NuTo::Brick8N::Brick8N] Exactly eight nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    mNodes[3] = rNodes[3];
    mNodes[4] = rNodes[4];
    mNodes[5] = rNodes[5];
    mNodes[6] = rNodes[6];
    mNodes[7] = rNodes[7];
    CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Brick8N::CalculateShapeFunctionsGeometry(const double rLocalCoordinates[3], std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==8);
    double plus_r,plus_s,plus_t,minus_r,minus_s,minus_t;

    plus_r  =1.0+rLocalCoordinates[0];
    plus_s  =1.0+rLocalCoordinates[1];
    plus_t  =1.0+rLocalCoordinates[2];

    minus_r =1.0-rLocalCoordinates[0];
    minus_s =1.0-rLocalCoordinates[1];
    minus_t =1.0-rLocalCoordinates[2];

    rShapeFunctions[0] =0.125 *minus_r *minus_s *minus_t;
    rShapeFunctions[1] =0.125 *plus_r  *minus_s *minus_t;
    rShapeFunctions[2] =0.125 *plus_r  *plus_s  *minus_t;
    rShapeFunctions[3] =0.125 *minus_r *plus_s  *minus_t;
    rShapeFunctions[4] =0.125 *minus_r *minus_s *plus_t;
    rShapeFunctions[5] =0.125 *plus_r  *minus_s *plus_t;
    rShapeFunctions[6] =0.125 *plus_r  *plus_s  *plus_t;
    rShapeFunctions[7] =0.125 *minus_r *plus_s  *plus_t;
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Brick8N::CalculateShapeFunctionsField(const double rLocalCoordinates[3], std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==8);
    double plus_r,plus_s,plus_t,minus_r,minus_s,minus_t;

    plus_r  =1.0+rLocalCoordinates[0];
    plus_s  =1.0+rLocalCoordinates[1];
    plus_t  =1.0+rLocalCoordinates[2];

    minus_r =1.0-rLocalCoordinates[0];
    minus_s =1.0-rLocalCoordinates[1];
    minus_t =1.0-rLocalCoordinates[2];

    rShapeFunctions[0] =0.125 *minus_r *minus_s *minus_t;
    rShapeFunctions[1] =0.125 *plus_r  *minus_s *minus_t;
    rShapeFunctions[2] =0.125 *plus_r  *plus_s  *minus_t;
    rShapeFunctions[3] =0.125 *minus_r *plus_s  *minus_t;
    rShapeFunctions[4] =0.125 *minus_r *minus_s *plus_t;
    rShapeFunctions[5] =0.125 *plus_r  *minus_s *plus_t;
    rShapeFunctions[6] =0.125 *plus_r  *plus_s  *plus_t;
    rShapeFunctions[7] =0.125 *minus_r *plus_s  *plus_t;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Brick8N::CalculateDerivativeShapeFunctionsGeometryNatural(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==24);
    double plus_r  =1.0+rLocalCoordinates[0];
    double plus_s  =1.0+rLocalCoordinates[1];
    double plus_t  =1.0+rLocalCoordinates[2];

    double minus_r =1.0-rLocalCoordinates[0];
    double minus_s =1.0-rLocalCoordinates[1];
    double minus_t =1.0-rLocalCoordinates[2];

    //node1
    rDerivativeShapeFunctions[0]  = -0.125 *minus_s *minus_t;
    rDerivativeShapeFunctions[1]  = -0.125 *minus_r *minus_t;
    rDerivativeShapeFunctions[2]  = -0.125 *minus_r *minus_s;
    //node2
    rDerivativeShapeFunctions[3]  =  0.125 *minus_s *minus_t;
    rDerivativeShapeFunctions[4]  = -0.125 *plus_r  *minus_t;
    rDerivativeShapeFunctions[5]  = -0.125 *plus_r  *minus_s;
    //node3
    rDerivativeShapeFunctions[6]  =  0.125 *plus_s  *minus_t;
    rDerivativeShapeFunctions[7]  =  0.125 *plus_r  *minus_t;
    rDerivativeShapeFunctions[8]  = -0.125 *plus_r  *plus_s;
    //node4
    rDerivativeShapeFunctions[9]  = -0.125 *plus_s  *minus_t;
    rDerivativeShapeFunctions[10] =  0.125 *minus_r *minus_t;
    rDerivativeShapeFunctions[11] = -0.125 *minus_r *plus_s;
    //node5
    rDerivativeShapeFunctions[12] = -0.125 *minus_s *plus_t;
    rDerivativeShapeFunctions[13] = -0.125 *minus_r *plus_t;
    rDerivativeShapeFunctions[14] =  0.125 *minus_r *minus_s;
    //node6
    rDerivativeShapeFunctions[15] =  0.125 *minus_s *plus_t;
    rDerivativeShapeFunctions[16] = -0.125 *plus_r  *plus_t;
    rDerivativeShapeFunctions[17] =  0.125 *plus_r  *minus_s;
    //node7
    rDerivativeShapeFunctions[18] =  0.125 *plus_s  *plus_t;
    rDerivativeShapeFunctions[19] =  0.125 *plus_r  *plus_t;
    rDerivativeShapeFunctions[20] =  0.125 *plus_r  *plus_s;
    //node8
    rDerivativeShapeFunctions[21] = -0.125 *plus_s  *plus_t;
    rDerivativeShapeFunctions[22] =  0.125 *minus_r *plus_t;
    rDerivativeShapeFunctions[23] =  0.125 *minus_r *plus_s;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Brick8N::CalculateDerivativeShapeFunctionsFieldNatural(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const
{
	CalculateDerivativeShapeFunctionsGeometryNatural(rLocalCoordinates,rDerivativeShapeFunctions);
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Brick8N::CalculateShapeFunctionsSurface(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==4);
    rShapeFunctions[0] = 0.25*(1.-rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
    rShapeFunctions[1] = 0.25*(1.+rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
    rShapeFunctions[2] = 0.25*(1.+rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
    rShapeFunctions[3] = 0.25*(1.-rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
}

//! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Brick8N::CalculateDerivativeShapeFunctionsLocalSurface(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
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

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Brick8N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{
	rSurfaceNodes.resize(4);
    std::cout << "check the surface node calculation for brick8 elements. It has not been verified." << std::endl;
	switch(rSurface)
    {
    case 0:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[1];
    	rSurfaceNodes[2] = mNodes[5];
    	rSurfaceNodes[3] = mNodes[4];
    	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[1];
    	rSurfaceNodes[1] = mNodes[2];
    	rSurfaceNodes[2] = mNodes[6];
    	rSurfaceNodes[3] = mNodes[5];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[3];
    	rSurfaceNodes[1] = mNodes[7];
    	rSurfaceNodes[2] = mNodes[6];
    	rSurfaceNodes[3] = mNodes[2];
    	break;
    case 3:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[4];
    	rSurfaceNodes[2] = mNodes[7];
    	rSurfaceNodes[3] = mNodes[3];
    	break;
    case 4:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[3];
    	rSurfaceNodes[2] = mNodes[2];
    	rSurfaceNodes[3] = mNodes[1];
    	break;
    case 5:
    	rSurfaceNodes[0] = mNodes[4];
    	rSurfaceNodes[1] = mNodes[5];
    	rSurfaceNodes[2] = mNodes[6];
    	rSurfaceNodes[3] = mNodes[7];
    	break;
    default:
    	throw MechanicsException("[NuTo::Tetrahedron4N::GetSurfaceNodes] there are only 4 surfaces for a tet4.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Brick8N::GetNumSurfaces()const
{
    return 4;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Brick8N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip;
}

//! @brief reorder element nodes
void NuTo::Brick8N::ReorderNodes()
{
    // swap 2 and 4
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[3];
    this->mNodes[3] = tmp;

    // swap 6 and 8
    tmp = this->mNodes[5];
    this->mNodes[5] = this->mNodes[7];
    this->mNodes[7] = tmp;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Brick8N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<8; count++)
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
template void NuTo::Brick8N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Brick8N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Brick8N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Brick8N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Brick8N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Brick8N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Brick8N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Brick8N" << std::endl;
#endif
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Solid)
           & BOOST_SERIALIZATION_NVP(mNodes);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Brick8N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Brick8N)
#endif // ENABLE_SERIALIZATION
