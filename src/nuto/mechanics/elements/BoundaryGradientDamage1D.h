// $Id: Truss.h 627 2013-05-22 07:43:22Z unger3 $
#ifndef BOUNDARYGRADIENTDAMAGE1D_H
#define BOUNDARYGRADIENTDAMAGE1D_H

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/Truss1D.h"

namespace NuTo
{
class StructureBase;
class DeformationGradient1D;
class ConstitutiveTangentLocal1x1;
class EngineeringStress1D;
class EngineeringStrain1D;
class HeatFlux1D;
class Damage;
class LocalEqPlasticStrain;
class NonlocalEqPlasticStrain;
template <int TNumRows, int TNumColumns> class ConstitutiveTangentLocal;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... boundary element for gradient models to ensure a more reasonable application of boundary conditions
class BoundaryGradientDamage1D : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    BoundaryGradientDamage1D(const StructureBase* rStructure,
    		std::vector<NuTo::NodeBase* >& rNodes,
    		Truss* rRealBoundaryElement,
    		bool rEdgeRealBoundaryElement,
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
        return NuTo::Element::BOUNDARYGRADIENTDAMAGE1D;
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes()const
    {
        return mNodes.size();
    }

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    int GetNumShapeFunctions()const
    {
    	return mNodes.size();
    }

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<(int)mNodes.size());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber)const
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<(int)mNodes.size());
        return mNodes[rLocalNodeNumber];
    }

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<(int)mNodes.size());
        mNodes[rLocalNodeNumber] = rNode;
    }

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
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
    {}

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    void GetIntegrationPointVolume(std::vector<double>& rVolume)const
    {}

    //! @brief returns the global coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    void GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
    {}

    //! @brief calculates the local coordinates of the nodes
    //! @param localCoordinates vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const;

    //! @brief stores the nonlocal eq plastic strain of the nodes
    //! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
    //! @param nonlocal eq plastic strain vector with already correct size allocated (2*nodes)
    //! this can be checked with an assertation
    void CalculateNodalNonlocalEqPlasticStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalEquivalentPlasticStrain)const;

    //! @brief stores the nonlocal total strain of the nodes
    //! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
    //! @param nonlocal total strain vector with already correct size allocated (1*nodes)
    //! this can be checked with an assertation
    void CalculateNodalNonlocalTotalStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalTotalStrain)const;

    //! @brief calculates output data fo the elmement
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    NuTo::Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput);

    //! @brief cast the base pointer, otherwise throws an exception
    const NuTo::BoundaryGradientDamage1D* AsBoundaryGradientDamage1D()const;

    //! @brief cast the base pointer, otherwise throws an exception
    NuTo::BoundaryGradientDamage1D* AsBoundaryGradientDamage1D();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief ... just for serialization
    BoundaryGradientDamage1D(){}

    //The real boundary element that is attached to the virtual boundary element
    const Truss* mRealBoundaryElement;


     std::vector<NodeBase*> mNodes;
    //edge of the real boundary element where the virtual boundary element is attached to
    //0 node 0
    //1 last node
    bool mEdgeRealBoundaryElement;

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    //! @param rNumDisp ... number of displacement dofs
    //! @param rNumTemp ... number of temperature dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs,int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqPlasticStrainDofs,int rNumNonlocalTotalStrainDofs) const;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    //! @param rNumDisp ... number of displacement dofs
    //! @param rNumTemp ... number of temperature dofs
    void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs,int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqPlasticStrainDofs,int rNumNonlocalTotalStrainDofs) const;

    //! @brief ... reorder nodes such that the sign of the length of the element changes
    void ReorderNodes()
    {
    }

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length is negative)
    void CheckElement()
    {
    }

};

} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::BoundaryGradientDamage1D)
#endif // ENABLE_SERIALIZATION

#endif //BOUNDARYGRADIENTDAMAGE1D_H
