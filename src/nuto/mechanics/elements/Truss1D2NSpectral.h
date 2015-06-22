#ifndef Truss1D2NSpectral_H
#define Truss1D2NSpectral_H


#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include <vector>
#include "nuto/mechanics/elements/Truss1D.h"
#include "nuto/mechanics/structures/StructureBase.h"

namespace NuTo
{

template <int TOrder>
class Truss1D2NSpectral : public Truss1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Truss1D2NSpectral(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
            ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType):
        NuTo::Truss1D::Truss1D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
    {
        if ((int)rNodes.size() != GetNumNodes())
            throw MechanicsException("[NuTo::Truss2NSpectral::Truss2NSpectral] Exactly (Order+1) nodes are required for this type of element.");

        for (int count=0; count<GetNumNodes(); count++)
            mNodes[count] = rNodes[count];

        this->CheckElement();

        //check the position of the nodes
        if (GetNumIntegrationPoints()!=GetNumNodes())
        {
            throw MechanicsException("[NuTo::Truss1D2NSpectral] Number of nodes and integration points must be identical for spectral elements.");
        }

        // TODO: Some calculations like in 2D??
    }

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
        switch (TOrder) {
        case 2:
            return NuTo::Element::TRUSS1D2NSPECTRALORDER2;
        case 3:
            return NuTo::Element::TRUSS1D2NSPECTRALORDER3;
        case 4:
            return NuTo::Element::TRUSS1D2NSPECTRALORDER4;
        default:
            throw MechanicsException("[NuTo::Truss1D2NSpectral] The specified order is not implemented yet.");
            break;
        }
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    inline int GetNumNodes()const
    {
        return TOrder+1;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    inline NodeBase* GetNode(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber) const
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < GetNumNodes());
        mNodes[rLocalNodeNumber] = rNode;
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesGeometry() const
    {
        return 2;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeGeometry(int rLocalNodeNumber) const
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < 2);
        switch(rLocalNodeNumber)
        {
            case 0:
                return mNodes[0];
                break;
            case 1:
                return mNodes[TOrder];
                break;
            default:
                throw MechanicsException("[Plane1D2NSpectralOrder2::GetNodeGeometry] element has only 2 nodes for geometry interpolation.");
        }
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeGeometry(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < 2);
        switch(rLocalNodeNumber)
        {
            case 0:
                return mNodes[0];
                break;
            case 1:
                return mNodes[TOrder];
                break;
            default:
                throw MechanicsException("[Plane1D2NSpectralOrder2::GetNodeGeometry] element has only 2 nodes for geometry interpolation.");
        }
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesField() const
    {
        return GetNumNodes();
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeField(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element (field interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeField(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
    {
        for (int count=0; count<GetNumNodes(); count++)
        {
            if (this->mNodes[count]==rOldPtr)
            {
                this->mNodes[count]=rNewPtr;
                break;
            }
        }
    }

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctionsNonlocalTotalStrain()const
    {
        return GetNumNodes();
    }

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
    {
        NuTo::ShapeFunctions1D::ShapeFunctionsTrussOrder1(rLocalCoordinates, rShapeFunctions);
    }

    //! @brief calculates the derivative of the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
    {
        NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1(rLocalCoordinates, rDerivativeShapeFunctions);
    }

    //! @brief calculates the shape functions for the field (2D is just the tensor product)
    //! @param rNaturalCoordinate natural coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsField(double rNaturalCoordinate, std::vector<double>& rShapeFunctions) const
    {
        switch (TOrder) {
        case 2:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussOrder2(rNaturalCoordinate, rShapeFunctions);
            break;
        case 3:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(rNaturalCoordinate, rShapeFunctions);
            break;
        case 4:
            NuTo::ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(rNaturalCoordinate, rShapeFunctions);
            break;
        default:
            break;
        }
    }

    //! @brief calculates the derivative of the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsField(double rNaturalCoordinate, std::vector<double>& rDerivativeShapeFunctions) const
    {
        switch (TOrder) {
        case 2:
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rNaturalCoordinate, rDerivativeShapeFunctions);
            break;
        case 3:
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder3(rNaturalCoordinate, rDerivativeShapeFunctions);
            break;
        case 4:
            NuTo::ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder4(rNaturalCoordinate, rDerivativeShapeFunctions);
            break;
        default:
            break;
        }
    }

    //! @brief returns the number of local degrees of freedom
    //! @return number of local degrees of freedom
    inline int GetNumLocalDofs() const
    {
        return GetNumNodes();
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeNonlocalTotalStrain(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeNonlocalTotalStrain(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber >= 0 && rLocalNodeNumber < GetNumNodes());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns the enum of the standard integration type for this element
    NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType()
    {
        switch (TOrder) {
        case 2:
            return NuTo::IntegrationType::IntegrationType1D2NLobatto3Ip;
        case 3:
            return NuTo::IntegrationType::IntegrationType1D2NLobatto4Ip;
        case 4:
            return NuTo::IntegrationType::IntegrationType1D2NLobatto5Ip;
        default:
            throw MechanicsException("[NuTo::Truss1D2NSpectral] The specified order is not implemented yet.");
            break;
        }
    }


protected:
    //! @brief ... just for serialization
    Truss1D2NSpectral(){}

    //! @brief ... reorder nodes such that the sign of the length of the element changes
    void ReorderNodes()
    {
        throw MechanicsException("[NuTo::Plane1D2NSpectral::ReorderNodes] not implemented.");
    }

    //! @brief element nodes
    NodeBase* mNodes[TOrder+1];
}; // class definition
} // namespace nuto

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane1D2NSpectral<2>)
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane1D2NSpectral<3>)
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane1D2NSpectral<4>)
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::Truss1D2NSpectral<2>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::Truss1D2NSpectral<3>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::Truss1D2NSpectral<4>)))
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif //TRUSS1D2NSPECTRAL_H
