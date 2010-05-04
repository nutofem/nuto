// $Id: $
#ifndef Truss1D3N_H
#define Truss1D3N_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <vector>
#include "nuto/mechanics/elements/Truss1D.h"

namespace NuTo
{

class Truss1D3N : public Truss1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Truss1D3N(NuTo::StructureBase* rStructure, std::vector<NuTo::NodeBase* >& rNodes,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Truss1D)
           & BOOST_SERIALIZATION_NVP(mNodes);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
        return NuTo::Element::TRUSS1D3N;
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes()const
    {
    	return 3;
    }

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctions()const
    {
    	return 3;
    }

    //! @brief calculates the shape functions
	//! @param rLocalCoordinates local coordinates of the integration point
	//! @param shape functions for all the nodes
	void CalculateShapeFunctions(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

	//! @brief calculates the derivative of the shape functions
	//! @param rLocalCoordinates local coordinates of the integration point
	//! @param derivative of the shape functions for all the nodes,
	//! first all the directions for a single node, and then for the next node
	void CalculateDerivativeShapeFunctions(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief returns the number of local degrees of freedom
    //! @return number of local degrees of freedom
    inline int GetNumLocalDofs()const
    {
    	return 3;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
    	assert(rLocalNodeNumber==0 || rLocalNodeNumber==1 || rLocalNodeNumber==2);
    	return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber)const
    {
    	assert(rLocalNodeNumber==0 || rLocalNodeNumber==1 || rLocalNodeNumber==2);
    	return mNodes[rLocalNodeNumber];
    }


    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
    	assert(rLocalNodeNumber==0 || rLocalNodeNumber==1 || rLocalNodeNumber==2);
    	mNodes[rLocalNodeNumber] = rNode;
    }

	//! @brief calculate list of global dofs related to the entries in the element stiffness matrix
    //! @param rGlobalDofsRow global dofs corresponding to the rows of the matrix
    //! @param rGlobalDofsColumn global dofs corresponding to the columns of the matrix
    void CalculateGlobalDofs(std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
    {
    	throw MechanicsException("[NuTo::Truss1D3N::CalculateGlobalDofs] to be implemented.");
    }

    //! @brief returns the enum of the standard integration type for this element
    NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType();


protected:
    //! @brief ... reorder nodes such that the sign of the length of the element changes
    void ReorderNodes();

    NodeBase* mNodes[3];
};
}
#endif //Truss1D3N_H
