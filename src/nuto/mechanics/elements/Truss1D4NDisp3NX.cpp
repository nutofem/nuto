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
#include "nuto/mechanics/elements/Truss1D4NDisp3NX.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>

NuTo::Truss1D4NDisp3NX::Truss1D4NDisp3NX(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Truss1D::Truss1D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=5)
    	throw MechanicsException("[NuTo::Truss1D4NDisp3NX::Truss1D4NDisp3NX] Exactly five nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
    mNodes[3] = rNodes[3];
    mNodes[4] = rNodes[4];
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D4NDisp3NX::CalculateShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    ShapeFunctions1D::ShapeFunctionsTrussOrder3(rLocalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D4NDisp3NX::CalculateShapeFunctionsNonlocalTotalStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
   ShapeFunctions1D::ShapeFunctionsTrussOrder2(rLocalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D4NDisp3NX::CalculateShapeFunctionsNonlocalEqStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    ShapeFunctions1D::ShapeFunctionsTrussOrder2(rLocalCoordinates, rShapeFunctions);
}


//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D4NDisp3NX::CalculateDerivativeShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder3(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D4NDisp3NX::CalculateDerivativeShapeFunctionsNonlocalTotalStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D4NDisp3NX::CalculateDerivativeShapeFunctionsNonlocalEqStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeGeometry(int rLocalNodeNumber)const
{
    assert(rLocalNodeNumber>=0 and rLocalNodeNumber<4 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[1];
    if (rLocalNodeNumber==2)
        return mNodes[3];
    if (rLocalNodeNumber==3)
        return mNodes[4];

    throw MechanicsException("[NuTo::Truss1D4NDisp3NX::GetNodeGeometry] the interpolation order is cubic for the geometry.");
    return 0;
}

//! @brief returns a pointer to the i-th node of the element (geometry interpolation)
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeGeometry(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber>=0 and rLocalNodeNumber<4 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[1];
    if (rLocalNodeNumber==2)
        return mNodes[3];
    if (rLocalNodeNumber==3)
        return mNodes[4];

    throw MechanicsException("[NuTo::Truss1D4NDisp3NX::GetNodeGeometry] the interpolation order is cubic for the geometry.");
    return 0;
}

//! @brief returns a pointer to the i-th node of the element (geometry interpolation)
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeField(int rLocalNodeNumber) const
{
    return GetNodeGeometry(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element (field interpolation)
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeField(int rLocalNodeNumber)
{
    return GetNodeGeometry(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeNonlocalTotalStrain(int rLocalNodeNumber)
{
    return GetNodeNonlocalEqStrain(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeNonlocalTotalStrain(int rLocalNodeNumber)const
{
    return GetNodeNonlocalEqStrain(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeNonlocalEqStrain(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber>=0 and rLocalNodeNumber<3 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[2];
    if (rLocalNodeNumber==2)
        return mNodes[4];

    throw MechanicsException("[NuTo::Truss1D4NDisp3NX::GetNodeNonlocalEqStrain] the interpolation order is quadratic for the nonlocal eq strain.");
    return 0;
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Truss1D4NDisp3NX::GetNodeNonlocalEqStrain(int rLocalNodeNumber)const
{
    assert(rLocalNodeNumber>=0 and rLocalNodeNumber<3 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[2];
    if (rLocalNodeNumber==2)
        return mNodes[4];

    throw MechanicsException("[NuTo::Truss1D4NDisp3NX::GetNodeNonlocalEqStrain] the interpolation order is quadratic for the nonlocal eq strain.");
    return 0;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Truss1D4NDisp3NX::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType1D2NGauss3Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Truss1D4NDisp3NX::ReorderNodes()
{
    std::cout << "reorder element nodes" << std::endl;
    // swap 0 and 4
    NodeBase* tmp = this->mNodes[0];
    this->mNodes[0] = this->mNodes[4];
    this->mNodes[4] = tmp;

    // swap 1 and 3
    tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[3];
    this->mNodes[3] = tmp;

}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Truss1D4NDisp3NX::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<5; count++)
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
template void NuTo::Truss1D4NDisp3NX::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D4NDisp3NX::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D4NDisp3NX::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D4NDisp3NX::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D4NDisp3NX::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D4NDisp3NX::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Truss1D4NDisp3NX::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Truss1D4NDisp3NX" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Truss1D)
       & BOOST_SERIALIZATION_NVP(mNodes);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Truss1D4NDisp3NX" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Truss1D4NDisp3NX)
#endif // ENABLE_SERIALIZATION
