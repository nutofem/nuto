// $Id: Plane2D3N.h 276 2010-06-30 13:04:32Z arnold2 $
#ifndef PLANE2D6N_H
#define PLANE2D6N_H

#include "nuto/mechanics/elements/Plane2D.h"

namespace NuTo
{

class Plane2D6N : public Plane2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Plane2D6N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes()const
    {
        return 6;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<6);
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber)const
    {
    	assert(rLocalNodeNumber>=0 && rLocalNodeNumber<6);
        return mNodes[rLocalNodeNumber];
    }


    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
    	assert(rLocalNodeNumber>=0 && rLocalNodeNumber<6);
        mNodes[rLocalNodeNumber] = rNode;
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesGeometry()const
    {
    	return 6;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeGeometry(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<6);
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeGeometry(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<6);
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesField()const
    {
    	return 6;
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeField(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<6);
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element (field interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeField(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<6);
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
        return NuTo::Element::PLANE2D6N;
    }

    //! @brief calculate the natural coordinates in 2D of all nodes
    static void CalculateNaturalNodeCoordinates(std::vector<std::array<double,2> >& rNaturalNodeCoordinates);

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr);

    //! @brief calculates the shape functions
    //! @param rNaturalCoordinates natural coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsGeometry(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the shape functions
    //! @param rNaturalCoordinates natural coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions
    //! @param rNaturalCoordinates natural coordinates (-1,1) of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions
    //! @param rNaturalCoordinates natural coordinates (-1,1) of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the shape functions for the surfaces (required for surface loads)
    //! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    void CalculateShapeFunctionsSurface(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const override;

    //! @brief calculates the derivative of the shape functions with respect to local coordinatesfor the surfaces (required for surface loads)
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsLocalSurface(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const override;

    //! @brief returns the surface nodes
    //! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
    //! @param surface nodes
    void GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const override;

    //! @brief returns the number of external surfaces
    //! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
    //! @param surface nodes
    int GetNumSurfaces()const override;

    //! @brief returns the enum of the standard integration type for this element
    NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType();


protected:
    //! @brief ... just for serialization
    Plane2D6N(){};
    //! @brief ... reorder nodes such that the sign of the length of the element changes
    void ReorderNodes();

    //! @brief element nodes
    NodeBase* mNodes[6];
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Plane2D6N)
#endif // ENABLE_SERIALIZATION

#endif //PLANE2D6N_H
