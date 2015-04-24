// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/Brick20N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include <assert.h>

//! @author Andrea Keszler, ISM
//! @date Octobre 2014
//! @brief ... brick element with 20 nodes - serendipitiy
NuTo::Brick20N::Brick20N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Solid::Solid(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
    if (rNodes.size()!=20)
        throw MechanicsException("[NuTo::Brick20N::Brick20N] Exactly twenty nodes are required for this type of element.");
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
    mNodes[15] = rNodes[15];
    mNodes[16] = rNodes[16];
    mNodes[17] = rNodes[17];
    mNodes[18] = rNodes[18];
    mNodes[19] = rNodes[19];
    mNodes[20] = rNodes[20];
    CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Brick20N::CalculateShapeFunctionsGeometry(const double rLocalCoordinates[3], std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==8);
    double plus_r,plus_s,plus_t,minus_r,minus_s,minus_t;

    plus_r  =1.0+rLocalCoordinates[0];
    plus_s  =1.0+rLocalCoordinates[1];
    plus_t  =1.0+rLocalCoordinates[2];

    minus_r =1.0-rLocalCoordinates[0];
    minus_s =1.0-rLocalCoordinates[1];
    minus_t =1.0-rLocalCoordinates[2];
    //for middle nodes $1-r^2$
    double mid_r =1.0-rLocalCoordinates[0]*rLocalCoordinates[0];
    double mid_s =1.0-rLocalCoordinates[1]*rLocalCoordinates[1];
    double mid_t =1.0-rLocalCoordinates[2]*rLocalCoordinates[2];

   // first 8 shape function are the same for Brick8N
    rShapeFunctions[0] =0.125 *minus_r *minus_s *minus_t;
    rShapeFunctions[1] =0.125 *plus_r  *minus_s *minus_t;
    rShapeFunctions[2] =0.125 *plus_r  *plus_s  *minus_t;
    rShapeFunctions[3] =0.125 *minus_r *plus_s  *minus_t;
    rShapeFunctions[4] =0.125 *minus_r *minus_s *plus_t;
    rShapeFunctions[5] =0.125 *plus_r  *minus_s *plus_t;
    rShapeFunctions[6] =0.125 *plus_r  *plus_s  *plus_t;
    rShapeFunctions[7] =0.125 *minus_r *plus_s  *plus_t;

    rShapeFunctions[8] =0.25*  mid_r*minus_s*minus_t;
    rShapeFunctions[9] =0.25* plus_r*  mid_s*minus_t;
    rShapeFunctions[10]=0.25*  mid_r* plus_s*minus_t;
    rShapeFunctions[11]=0.25*minus_r*  mid_s*minus_t;

    rShapeFunctions[12]=0.25*minus_r*minus_s*  mid_t;
    rShapeFunctions[13]=0.25* plus_r*minus_s*  mid_t;
    rShapeFunctions[14]=0.25* plus_r* plus_s*  mid_t;
    rShapeFunctions[15]=0.25*minus_r* plus_s*  mid_t;

    rShapeFunctions[16]=0.25*  mid_r*minus_s* plus_t;
    rShapeFunctions[17]=0.25* plus_r*  mid_s* plus_t;
    rShapeFunctions[18]=0.25*  mid_r* plus_s* plus_t;
    rShapeFunctions[19]=0.25*minus_r*  mid_s* plus_t;
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Brick20N::CalculateShapeFunctionsField(const double rLocalCoordinates[3], std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==8);
    double plus_r,plus_s,plus_t,minus_r,minus_s,minus_t;

    plus_r  =1.0+rLocalCoordinates[0];
    plus_s  =1.0+rLocalCoordinates[1];
    plus_t  =1.0+rLocalCoordinates[2];

    minus_r =1.0-rLocalCoordinates[0];
    minus_s =1.0-rLocalCoordinates[1];
    minus_t =1.0-rLocalCoordinates[2];
    //for middle nodes $1-r^2$
    double mid_r =1.0-rLocalCoordinates[0]*rLocalCoordinates[0];
    double mid_s =1.0-rLocalCoordinates[1]*rLocalCoordinates[1];
    double mid_t =1.0-rLocalCoordinates[2]*rLocalCoordinates[2];

   // first 8 shape function are the same for Brick8N
    rShapeFunctions[0] =0.125 *minus_r *minus_s *minus_t;
    rShapeFunctions[1] =0.125 *plus_r  *minus_s *minus_t;
    rShapeFunctions[2] =0.125 *plus_r  *plus_s  *minus_t;
    rShapeFunctions[3] =0.125 *minus_r *plus_s  *minus_t;
    rShapeFunctions[4] =0.125 *minus_r *minus_s *plus_t;
    rShapeFunctions[5] =0.125 *plus_r  *minus_s *plus_t;
    rShapeFunctions[6] =0.125 *plus_r  *plus_s  *plus_t;
    rShapeFunctions[7] =0.125 *minus_r *plus_s  *plus_t;

    rShapeFunctions[8] =0.25*  mid_r*minus_s*minus_t;
    rShapeFunctions[9] =0.25* plus_r*  mid_s*minus_t;
    rShapeFunctions[10]=0.25*  mid_r* plus_s*minus_t;
    rShapeFunctions[11]=0.25*minus_r*  mid_s*minus_t;

    rShapeFunctions[12]=0.25*minus_r*minus_s*  mid_t;
    rShapeFunctions[13]=0.25* plus_r*minus_s*  mid_t;
    rShapeFunctions[14]=0.25* plus_r* plus_s*  mid_t;
    rShapeFunctions[15]=0.25*minus_r* plus_s*  mid_t;

    rShapeFunctions[16]=0.25*  mid_r*minus_s* plus_t;
    rShapeFunctions[17]=0.25* plus_r*  mid_s* plus_t;
    rShapeFunctions[18]=0.25*  mid_r* plus_s* plus_t;
    rShapeFunctions[19]=0.25*minus_r*  mid_s* plus_t;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Brick20N::CalculateDerivativeShapeFunctionsGeometryNatural(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==60);
    double plus_r  =1.0+rLocalCoordinates[0];
    double plus_s  =1.0+rLocalCoordinates[1];
    double plus_t  =1.0+rLocalCoordinates[2];

    double minus_r =1.0-rLocalCoordinates[0];
    double minus_s =1.0-rLocalCoordinates[1];
    double minus_t =1.0-rLocalCoordinates[2];

    //for middle nodes $1-r^2$
    double mid_r =1.0-rLocalCoordinates[0]*rLocalCoordinates[0];
    double mid_s =1.0-rLocalCoordinates[1]*rLocalCoordinates[1];
    double mid_t =1.0-rLocalCoordinates[2]*rLocalCoordinates[2];
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
    //node9
    rDerivativeShapeFunctions[24] =+0.25*minus_s*minus_t;
    rDerivativeShapeFunctions[25] =-0.25*  mid_r*minus_t;
    rDerivativeShapeFunctions[26] =-0.25*  mid_r*minus_s;
    //node10
    rDerivativeShapeFunctions[27] =+0.25*  mid_s*minus_t;
    rDerivativeShapeFunctions[28] =+0.25* plus_r*minus_t;
    rDerivativeShapeFunctions[29] =-0.25* plus_r*  mid_s;
    //node11
    rDerivativeShapeFunctions[30]=+0.25* plus_s*minus_t;
    rDerivativeShapeFunctions[31]=+0.25*  mid_r*minus_t;
    rDerivativeShapeFunctions[32]=-0.25*  mid_r* plus_s;
    //node12
    rDerivativeShapeFunctions[33]=-0.25*  mid_s*minus_t;
    rDerivativeShapeFunctions[34]=+0.25*minus_r*minus_t;
    rDerivativeShapeFunctions[35]=-0.25*minus_r*  mid_s;
    //node13
    rDerivativeShapeFunctions[36]=0.25*minus_s*  mid_t;
    rDerivativeShapeFunctions[37]=0.25*minus_r*  mid_t;
    rDerivativeShapeFunctions[38]=0.25*minus_r*minus_s;
    //node14
    rDerivativeShapeFunctions[39]=+0.25*minus_s*  mid_t;
    rDerivativeShapeFunctions[40]=-0.25* plus_r*  mid_t;
    rDerivativeShapeFunctions[41]=+0.25* plus_r*minus_s;
    //node15
    rDerivativeShapeFunctions[42]=+0.25* plus_s*  mid_t;
    rDerivativeShapeFunctions[43]=+0.25* plus_r*  mid_t;
    rDerivativeShapeFunctions[44]=+0.25* plus_r* plus_s;
    //node16
    rDerivativeShapeFunctions[45]=-0.25* plus_s*  mid_t;
    rDerivativeShapeFunctions[46]=+0.25*minus_r*  mid_t;
    rDerivativeShapeFunctions[47]=+0.25*minus_r* plus_s;
    //node17
    rDerivativeShapeFunctions[48]=+0.25*minus_s* plus_t;
    rDerivativeShapeFunctions[49]=-0.25*  mid_r* plus_t;
    rDerivativeShapeFunctions[50]=+0.25*  mid_r*minus_s;
    //node18
    rDerivativeShapeFunctions[51]=+0.25*  mid_s* plus_t;
    rDerivativeShapeFunctions[52]=+0.25* plus_r* plus_t;
    rDerivativeShapeFunctions[53]=+0.25* plus_r*  mid_s;
    //node19
    rDerivativeShapeFunctions[54]=+0.25* plus_s* plus_t;
    rDerivativeShapeFunctions[55]=+0.25*  mid_r* plus_t;
    rDerivativeShapeFunctions[56]=+0.25*  mid_r* plus_s;
    //node20
    rDerivativeShapeFunctions[57]=-0.25*  mid_s* plus_t;
    rDerivativeShapeFunctions[58]=+0.25*minus_r* plus_t;
    rDerivativeShapeFunctions[59]=+0.25*minus_r*  mid_s;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Brick20N::CalculateDerivativeShapeFunctionsFieldNatural(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const
{
	CalculateDerivativeShapeFunctionsGeometryNatural(rLocalCoordinates,rDerivativeShapeFunctions);
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Brick20N::CalculateShapeFunctionsSurface(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==8);
	// surface 2d definitions!!!
	double minus_r=1.-rNaturalCoordinates[0];
	double  plus_r=1.+rNaturalCoordinates[0];
	double minus_s=1.-rNaturalCoordinates[1];
	double  plus_s=1.+rNaturalCoordinates[1];

	rShapeFunctions[0] = 0.25*minus_r*minus_s;
    rShapeFunctions[1] = 0.25* plus_r*minus_s;
    rShapeFunctions[2] = 0.25* plus_r* plus_s;
    rShapeFunctions[3] = 0.25*minus_r* plus_s;

    //for middle nodes $1-r^2$
    double mid_r =1.0-rNaturalCoordinates[0]*rNaturalCoordinates[0];
    double mid_s =1.0-rNaturalCoordinates[1]*rNaturalCoordinates[1];
    rShapeFunctions[4] =0.25*mid_r*minus_s;
    rShapeFunctions[5] =0.25* plus_r*  mid_s;
    rShapeFunctions[6] =0.25*  mid_r* plus_s;
    rShapeFunctions[7] =0.25*minus_r*  mid_s;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates for the surfaces (required for surface loads)
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Brick20N::CalculateDerivativeShapeFunctionsLocalSurface(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==16);
    rDerivativeShapeFunctions[0] = -0.25*(1.-rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[1] = -0.25*(1.-rNaturalCoordinates[0]);

    rDerivativeShapeFunctions[2] = +0.25*(1.-rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[3] = -0.25*(1.+rNaturalCoordinates[0]);

    rDerivativeShapeFunctions[4] = +0.25*(1.+rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[5] = +0.25*(1.+rNaturalCoordinates[0]);

    rDerivativeShapeFunctions[6] = -0.25*(1.+rNaturalCoordinates[1]);
    rDerivativeShapeFunctions[7] = +0.25*(1.-rNaturalCoordinates[0]);

	double minus_r=1.-rNaturalCoordinates[0];
	double  plus_r=1.+rNaturalCoordinates[0];
	double minus_s=1.-rNaturalCoordinates[1];
	double  plus_s=1.+rNaturalCoordinates[1];

	double mid_r =1.0-rNaturalCoordinates[0]*rNaturalCoordinates[0];
    double mid_s =1.0-rNaturalCoordinates[1]*rNaturalCoordinates[1];

	rDerivativeShapeFunctions[0] = -0.25*minus_s;
	rDerivativeShapeFunctions[1] = -0.25*minus_r;

    rDerivativeShapeFunctions[2] = +0.25* minus_s;
    rDerivativeShapeFunctions[3] = -0.25* plus_r;

    rDerivativeShapeFunctions[4] = +0.25* plus_s;
    rDerivativeShapeFunctions[5] = +0.25* plus_r;

    rDerivativeShapeFunctions[6] = -0.25* plus_s;
    rDerivativeShapeFunctions[7] = +0.25*minus_r;

    rDerivativeShapeFunctions[8] =+0.25*minus_s;
    rDerivativeShapeFunctions[9] =-0.25*mid_r;

    rDerivativeShapeFunctions[10] =+0.25* mid_s;
    rDerivativeShapeFunctions[11] =+0.25* plus_r;

    rDerivativeShapeFunctions[12] =+0.25* plus_s;
    rDerivativeShapeFunctions[13] =+0.25*  mid_r;

    rDerivativeShapeFunctions[14] =-0.25*  mid_s;
    rDerivativeShapeFunctions[15] =+0.25*minus_r;

}

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Brick20N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
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

      	rSurfaceNodes[4] = mNodes[8];
       	rSurfaceNodes[5] = mNodes[13];
       	rSurfaceNodes[6] = mNodes[16];
       	rSurfaceNodes[7] = mNodes[12];
   	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[1];
    	rSurfaceNodes[1] = mNodes[2];
    	rSurfaceNodes[2] = mNodes[6];
    	rSurfaceNodes[3] = mNodes[5];

      	rSurfaceNodes[4] = mNodes[9];
       	rSurfaceNodes[5] = mNodes[14];
       	rSurfaceNodes[6] = mNodes[17];
       	rSurfaceNodes[7] = mNodes[13];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[3];
    	rSurfaceNodes[1] = mNodes[7];
    	rSurfaceNodes[2] = mNodes[6];
    	rSurfaceNodes[3] = mNodes[2];

      	rSurfaceNodes[4] = mNodes[15];
       	rSurfaceNodes[5] = mNodes[18];
       	rSurfaceNodes[6] = mNodes[14];
       	rSurfaceNodes[7] = mNodes[10];
   	break;
    case 3:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[4];
    	rSurfaceNodes[2] = mNodes[7];
    	rSurfaceNodes[3] = mNodes[3];

      	rSurfaceNodes[4] = mNodes[13];
       	rSurfaceNodes[5] = mNodes[19];
       	rSurfaceNodes[6] = mNodes[15];
       	rSurfaceNodes[7] = mNodes[11];
    	break;
    case 4:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[3];
    	rSurfaceNodes[2] = mNodes[2];
    	rSurfaceNodes[3] = mNodes[1];

       	rSurfaceNodes[4] = mNodes[11];
       	rSurfaceNodes[5] = mNodes[5];
       	rSurfaceNodes[6] = mNodes[9];
       	rSurfaceNodes[7] = mNodes[8];
    	break;
    case 5:
    	rSurfaceNodes[0] = mNodes[4];
    	rSurfaceNodes[1] = mNodes[5];
    	rSurfaceNodes[2] = mNodes[6];
    	rSurfaceNodes[3] = mNodes[7];

       	rSurfaceNodes[4] = mNodes[16];
       	rSurfaceNodes[5] = mNodes[17];
       	rSurfaceNodes[6] = mNodes[18];
       	rSurfaceNodes[7] = mNodes[19];
     	break;
    default:
    	throw MechanicsException("[NuTo::Brick20N::GetSurfaceNodes] there are only 6 surfaces for a brick element.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
inline int NuTo::Brick20N::GetNumSurfaces()const
{
    return 6;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Brick20N::GetStandardIntegrationType()
{
    //! @todo add Integration type IntegrationType3D20NGauss3x3x3Ip for Brick 20N
	std::cout<<" Attention: Brick20N needs Gauss 3x3x3 integration, here 2x2x2 so far.\n";
	return NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip;
//	return NuTo::IntegrationType::IntegrationType3D20NGauss3x3x3Ip;
}

//! @brief reorder element nodes
void NuTo::Brick20N::ReorderNodes()
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
void NuTo::Brick20N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<20; count++)
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
template void NuTo::Brick20N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Brick20N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Brick20N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Brick20N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Brick20N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Brick20N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Brick20N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Brick20N" << std::endl;
#endif
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Solid)
           & BOOST_SERIALIZATION_NVP(mNodes);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Brick20N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Brick20N)
#endif // ENABLE_SERIALIZATION
