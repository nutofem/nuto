// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/Tetrahedron10N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>

NuTo::Tetrahedron10N::Tetrahedron10N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Solid::Solid(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
    if (rNodes.size()!=10)
        throw MechanicsException("[NuTo::Tetrahedron10N::Tetrahedron10N] Exactly ten nodes are required for this type of element.");
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
    CheckElement();
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Tetrahedron10N::CalculateShapeFunctionsGeometry(const double rLocalCoordinates[3], std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==10);
    double r,s,t;

    r  =rLocalCoordinates[0];
    s  =rLocalCoordinates[1];
    t  =rLocalCoordinates[2];

    // node 0 (0,0,0)
    rShapeFunctions[0] =1.-3.*(r+s+t)+2.*(r*r+s*s+t*t)+4.*(r*s+r*t+s*t);

    // node 1 (1,0,0)
    rShapeFunctions[1] =-r+2.*r*r;

    // node 2 (0,1,0)
    rShapeFunctions[2] =-s+2.*s*s;

    // node 3 (0,0,1)
    rShapeFunctions[3] =-t+2.*t*t;

    // node 4 (0.5,0,0)
    rShapeFunctions[4] =4.*r*(1.-r-s-t);

    // node 5 (0.5,0.5,0)
    rShapeFunctions[5] =4.*r*s;

    // node 6 (0,0.5,0)
    rShapeFunctions[6] =4.*s*(1.-r-s-t);

    // node 7 (0,0,0.5)
    rShapeFunctions[7] =4.*t*(1.-r-s-t);

    // node 8 (0,0.5,0.5)
    rShapeFunctions[8] =4.*s*t;

    // node 9 (0.5,0,0.5)
    rShapeFunctions[9] =4.*r*t;
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Tetrahedron10N::CalculateShapeFunctionsField(const double rLocalCoordinates[3], std::vector<double>& rShapeFunctions)const
{
    assert(rShapeFunctions.size()==10);
    double r,s,t;

    r  =rLocalCoordinates[0];
    s  =rLocalCoordinates[1];
    t  =rLocalCoordinates[2];

    // node 0 (0,0,0)
    rShapeFunctions[0] =1.-3.*(r+s+t)+2.*(r*r+s*s+t*t)+4.*(r*s+r*t+s*t);

    // node 1 (1,0,0)
    rShapeFunctions[1] =-r+2.*r*r;

    // node 2 (0,1,0)
    rShapeFunctions[2] =-s+2.*s*s;

    // node 3 (0,0,1)
    rShapeFunctions[3] =-t+2.*t*t;

    // node 4 (0.5,0,0)
    rShapeFunctions[4] =4.*r*(1.-r-s-t);

    // node 5 (0.5,0.5,0)
    rShapeFunctions[5] =4.*r*s;

    // node 6 (0,0.5,0)
    rShapeFunctions[6] =4.*s*(1.-r-s-t);

    // node 7 (0,0,0.5)
    rShapeFunctions[7] =4.*t*(1.-r-s-t);

    // node 8 (0,0.5,0.5)
    rShapeFunctions[8] =4.*s*t;

    // node 9 (0.5,0,0.5)
    rShapeFunctions[9] =4.*r*t;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Tetrahedron10N::CalculateDerivativeShapeFunctionsGeometryNatural(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==30);
    double r,s,t;

    r  =rLocalCoordinates[0];
    s  =rLocalCoordinates[1];
    t  =rLocalCoordinates[2];

    //node1
    rDerivativeShapeFunctions[0]  = -3.+ 4.*(r+s+t);
    rDerivativeShapeFunctions[1]  = -3.+ 4.*(r+s+t);
    rDerivativeShapeFunctions[2]  = -3.+ 4.*(r+s+t);

    //node2
    rDerivativeShapeFunctions[3]  = -1.+4.*r;
    rDerivativeShapeFunctions[4]  =  0;
    rDerivativeShapeFunctions[5]  =  0;

    //node3
    rDerivativeShapeFunctions[6]  =  0;
    rDerivativeShapeFunctions[7]  = -1.+4.*s;
    rDerivativeShapeFunctions[8]  =  0;

    //node4
    rDerivativeShapeFunctions[9]  =  0;
    rDerivativeShapeFunctions[10] =  0;
    rDerivativeShapeFunctions[11] = -1.+4.*t ;

    //node5
    rDerivativeShapeFunctions[12] =  4.-8.*r-4.*s-4.*t;
    rDerivativeShapeFunctions[13] = -4.*r;
    rDerivativeShapeFunctions[14] = -4.*r;

    //node6
    rDerivativeShapeFunctions[15] =  4.*s;
    rDerivativeShapeFunctions[16] =  4.*r;
    rDerivativeShapeFunctions[17] =  0.;

    //node7
    rDerivativeShapeFunctions[18] = -4.*s ;
    rDerivativeShapeFunctions[19] =  4.-8.*s-4.*r-4.*t;
    rDerivativeShapeFunctions[20] = -4.*s;

    //node8
    rDerivativeShapeFunctions[21] = -4.*t;
    rDerivativeShapeFunctions[22] = -4.*t;
    rDerivativeShapeFunctions[23] =  4.-8.*t-4.*r-4.*s;

    //node9
    rDerivativeShapeFunctions[24] =  0.;
    rDerivativeShapeFunctions[25] =  4.*t;
    rDerivativeShapeFunctions[26] =  4.*s;

    //node10
    rDerivativeShapeFunctions[27] =  4.*t;
    rDerivativeShapeFunctions[28] =  0;
    rDerivativeShapeFunctions[29] =  4.*r;
}

//! @brief calculates the derivative of the shape functions with respect to local coordinates
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Tetrahedron10N::CalculateDerivativeShapeFunctionsFieldNatural(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const
{
    assert(rDerivativeShapeFunctions.size()==30);
    double r,s,t;

    r  =rLocalCoordinates[0];
    s  =rLocalCoordinates[1];
    t  =rLocalCoordinates[2];

    //node1
    rDerivativeShapeFunctions[0]  = -3.+ 4.*(r+s+t);
    rDerivativeShapeFunctions[1]  = -3.+ 4.*(r+s+t);
    rDerivativeShapeFunctions[2]  = -3.+ 4.*(r+s+t);

    //node2
    rDerivativeShapeFunctions[3]  = -1.+4.*r;
    rDerivativeShapeFunctions[4]  =  0;
    rDerivativeShapeFunctions[5]  =  0;

    //node3
    rDerivativeShapeFunctions[6]  =  0;
    rDerivativeShapeFunctions[7]  = -1.+4.*s;
    rDerivativeShapeFunctions[8]  =  0;

    //node4
    rDerivativeShapeFunctions[9]  =  0;
    rDerivativeShapeFunctions[10] =  0;
    rDerivativeShapeFunctions[11] = -1.+4.*t ;

    //node5
    rDerivativeShapeFunctions[12] =  4.-8.*r-4.*s-4.*t;
    rDerivativeShapeFunctions[13] = -4.*r;
    rDerivativeShapeFunctions[14] = -4.*r;

    //node6
    rDerivativeShapeFunctions[15] =  4.*s;
    rDerivativeShapeFunctions[16] =  4.*r;
    rDerivativeShapeFunctions[17] =  0.;

    //node7
    rDerivativeShapeFunctions[18] = -4.*s ;
    rDerivativeShapeFunctions[19] =  4.-8.*s-4.*r-4.*t;
    rDerivativeShapeFunctions[20] = -4.*s;

    //node8
    rDerivativeShapeFunctions[21] = -4.*t;
    rDerivativeShapeFunctions[22] = -4.*t;
    rDerivativeShapeFunctions[23] =  4.-8.*t-4.*r-4.*s;

    //node9
    rDerivativeShapeFunctions[24] =  0.;
    rDerivativeShapeFunctions[25] =  4.*t;
    rDerivativeShapeFunctions[26] =  4.*s;

    //node10
    rDerivativeShapeFunctions[27] =  4.*t;
    rDerivativeShapeFunctions[28] =  0;
    rDerivativeShapeFunctions[29] =  4.*r;
}

//! @brief calculates the shape functions for the surfaces (required for surface loads)
//! @param rNaturalCoordinates local coordinates of the integration point (in the local surface coordinate system)
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Tetrahedron10N::CalculateShapeFunctionsSurface(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
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

//! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
//! @param rNaturalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Tetrahedron10N::CalculateDerivativeShapeFunctionsLocalSurface(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
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

    rDerivativeShapeFunctions[10] = -4.*s;;
    rDerivativeShapeFunctions[11] = -8.*s-4.*r+4.;


}

//! @brief returns the surface nodes
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
void NuTo::Tetrahedron10N::GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const
{
	rSurfaceNodes.resize(6);
	switch(rSurface)
    {
    case 0:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[1];
    	rSurfaceNodes[2] = mNodes[3];
    	rSurfaceNodes[3] = mNodes[4];
    	rSurfaceNodes[4] = mNodes[9];
    	rSurfaceNodes[5] = mNodes[7];
    	break;
    case 1:
    	rSurfaceNodes[0] = mNodes[2];
    	rSurfaceNodes[1] = mNodes[3];
    	rSurfaceNodes[2] = mNodes[1];
    	rSurfaceNodes[3] = mNodes[8];
    	rSurfaceNodes[4] = mNodes[9];
    	rSurfaceNodes[5] = mNodes[5];
    	break;
    case 2:
    	rSurfaceNodes[0] = mNodes[1];
    	rSurfaceNodes[1] = mNodes[0];
    	rSurfaceNodes[2] = mNodes[2];
    	rSurfaceNodes[3] = mNodes[4];
    	rSurfaceNodes[4] = mNodes[6];
    	rSurfaceNodes[5] = mNodes[5];
    	break;
    case 3:
    	rSurfaceNodes[0] = mNodes[0];
    	rSurfaceNodes[1] = mNodes[3];
    	rSurfaceNodes[2] = mNodes[2];
    	rSurfaceNodes[3] = mNodes[7];
    	rSurfaceNodes[4] = mNodes[8];
    	rSurfaceNodes[5] = mNodes[6];
    	break;
    default:
    	throw MechanicsException("[NuTo::Tetrahedron4N::GetSurfaceNodes] there are only 4 surfaces for a tet4.");
    }

}

//! @brief returns the number of external surfaces
//! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
//! @param surface nodes
int NuTo::Tetrahedron10N::GetNumSurfaces()const
{
    return 4;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Tetrahedron10N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType3D4NGauss4Ip;
}


//! @brief reorder element nodes
void NuTo::Tetrahedron10N::ReorderNodes()
{
    // swap nodes 1 and 2
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[2];
    this->mNodes[2] = tmp;

    // swap nodes 4 and 6
    tmp = this->mNodes[4];
    this->mNodes[4] = this->mNodes[6];
    this->mNodes[6] = tmp;

    // swap nodes 8 and 9
    tmp = this->mNodes[8];
    this->mNodes[8] = this->mNodes[9];
    this->mNodes[9] = tmp;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Tetrahedron10N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
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
template void NuTo::Tetrahedron10N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron10N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron10N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron10N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron10N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Tetrahedron10N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Tetrahedron10N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Tetrahedron10N" << std::endl;
#endif
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Solid)
           & BOOST_SERIALIZATION_NVP(mNodes);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Tetrahedron10N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Tetrahedron10N)
#endif // ENABLE_SERIALIZATION
