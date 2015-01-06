// $Id$
#ifndef TRUSS1D4NDISP3NX_H
#define TRUSS1D4NDISP3NX_H

#include "nuto/mechanics/elements/Truss1D.h"

namespace NuTo
{

class Truss1D4NDisp3NX : public Truss1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Truss1D4NDisp3NX(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
        return NuTo::Element::TRUSS1D4NDISP3NX;
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes()const
    {
        return 5;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<5);
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber)const
    {
    	assert(rLocalNodeNumber>=0 && rLocalNodeNumber<5);
        return mNodes[rLocalNodeNumber];
    }


    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
    	assert(rLocalNodeNumber>=0 && rLocalNodeNumber<5);
        mNodes[rLocalNodeNumber] = rNode;
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesGeometry()const
    {
    	return 4;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeGeometry(int rLocalNodeNumber)const;

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeGeometry(int rLocalNodeNumber);

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesField()const
    {
    	return 4;
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeField(int rLocalNodeNumber) const;

    //! @brief returns a pointer to the i-th node of the element (field interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeField(int rLocalNodeNumber);

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr);

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    int GetNumShapeFunctionsNonlocalTotalStrain()const
    {
    	return 3;
    }

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    int GetNumShapeFunctionsNonlocalEqStrain()const
    {
        return 3;
    }

    //! @brief calculates the shape functions
	//! @param rLocalCoordinates local coordinates of the integration point
	//! @param shape functions for all the nodes
	void CalculateShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsNonlocalTotalStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsNonlocalEqStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

	//! @brief calculates the derivative of the shape functions
	//! @param rLocalCoordinates local coordinates of the integration point
	//! @param derivative of the shape functions for all the nodes,
	//! first all the directions for a single node, and then for the next node
	void CalculateDerivativeShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsNonlocalTotalStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsNonlocalEqStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeNonlocalTotalStrain(int rLocalNodeNumber);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeNonlocalTotalStrain(int rLocalNodeNumber)const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeNonlocalEqStrain(int rLocalNodeNumber);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeNonlocalEqStrain(int rLocalNodeNumber)const;

	//! @brief calculate list of global dofs related to the entries in the element stiffness matrix
    //! @param rGlobalDofsRow global dofs corresponding to the rows of the matrix
    //! @param rGlobalDofsColumn global dofs corresponding to the columns of the matrix
    void CalculateGlobalDofs(std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
    {
    	throw MechanicsException("[NuTo::Truss1D4NDisp3NX::CalculateGlobalDofs] to be implemented.");
    }

    //! @brief returns the enum of the standard integration type for this element
    NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType();


protected:
    //! @brief ... just for serialization
    Truss1D4NDisp3NX(){}

    //! @brief ... reorder nodes such that the sign of the length of the element changes
    void ReorderNodes();

    NodeBase* mNodes[5];
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Truss1D4NDisp3NX)
#endif // ENABLE_SERIALIZATION
#endif //TRUSS1D4NDISP3NX_H
