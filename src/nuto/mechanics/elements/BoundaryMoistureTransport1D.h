#ifndef BOUNDARYMOISTURETRANSPORT1D_H
#define BOUNDARYMOISTURETRANSPORT1D_H

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/Truss1D.h"

#include "nuto/mechanics/elements/BoundaryGradientDamage1D.h"

namespace NuTo
{

//! @author Volker Hirthammer, BAM
//! @date March 2015
//! @brief ... boundary element for moisture transport model
class BoundaryMoistureTransport1D : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //@brief: constructor

    BoundaryMoistureTransport1D(const StructureBase* rStructure,
            Truss* rRealBoundaryElement,
            int rSurfaceEdge,
            BoundaryCondition::eType rBoundaryConditionType,
            ElementData::eElementDataType rElementDataType,
            IntegrationType::eIntegrationType rIntegrationType,
            IpData::eIpDataType rIpDataType
            );


    BoundaryMoistureTransport1D(const StructureBase* rStructure,
            Truss* rRealBoundaryElement,
            NodeBase* rSurfaceNode,
            BoundaryCondition::eType rBoundaryConditionType,
            ElementData::eElementDataType rElementDataType,
            IntegrationType::eIntegrationType rIntegrationType,
            IpData::eIpDataType rIpDataType
            );

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension()const
    {
        return 1;
    }

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetLocalDimension()const
    {
        return 1;
    }

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
        return NuTo::Element::BOUNDARYMOISTURETRANSPORT1D;
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes()const
    {
        return 1;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber==0);
        std::vector<const NodeBase*> surfaceNodes(1);
        mRealBoundaryElement->GetSurfaceNodes(mSurfaceEdge, surfaceNodes);
        return const_cast<NodeBase*>(surfaceNodes[0]);
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber==0);
        std::vector<const NodeBase*> surfaceNodes(1);
        mRealBoundaryElement->GetSurfaceNodes(mSurfaceEdge, surfaceNodes);
        return surfaceNodes[0];
    }


    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
        // todo?
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesGeometry()const
    {
        return 1;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeGeometry(int rLocalNodeNumber)const
    {
        return GetNode(rLocalNodeNumber);
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeGeometry(int rLocalNodeNumber)
    {
        return GetNode(rLocalNodeNumber);
    }

    //! @brief returns the number of nodes in this element(geometry interpolation)
    //! @return number of nodes
    int GetNumNodesField()const
    {
        return 1;
    }

    //! @brief returns a pointer to the i-th node of the element (geometry interpolation)
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNodeField(int rLocalNodeNumber)const
    {
        return GetNode(rLocalNodeNumber);
    }

    //! @brief returns a pointer to the i-th node of the element (field interpolation)
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNodeField(int rLocalNodeNumber)
    {
        return GetNode(rLocalNodeNumber);
    }

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr);

    //! @brief sets the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @param rSection pointer to section
    //! @return pointer to constitutive law
    void SetSection(const SectionBase* rSection)
    {
        throw MechanicsException("[NuTo::BoundaryGradientDamage1D::SetSection] The section is defined via the real boundary element.");
    }

    //! @brief returns a pointer to the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @return pointer to section
    const SectionBase* GetSection()const
    {
        return mRealBoundaryElement->GetSection();
    }

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    void GetIntegrationPointVolume(std::vector<double>& rVolume)const
    {}

    //! @brief returns the global coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    void GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
    {}

    //! @brief calculates output data fo the elmement
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    NuTo::Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput);

protected:
    //! @brief ... just for serialization
    BoundaryMoistureTransport1D(){}

    //The real boundary element that is attached to the virtual boundary element
    const Truss* mRealBoundaryElement;

    // edge number 0.. left, 1.. right
    int mSurfaceEdge;

    BoundaryCondition::eType mBoundaryConditionType;


    // build global row dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int rNumRelativeHumidityDofs, int rNumWaterPhaseFractionDofs) const;


    //! @brief returns the mSurfaceEdge variable for the node rNode
    int CalculateSurfaceEdge(const NodeBase* rNode) const;

    //! @brief ... reorder nodes such that the sign of the length of the element changes
    void ReorderNodes()
    {}

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length is negative)
    void CheckElement()
    {}

};

} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::BoundaryGradientDamage1D)
#endif // ENABLE_SERIALIZATION


#endif // BOUNDARYMOISTURETRANSPORT1D_H
