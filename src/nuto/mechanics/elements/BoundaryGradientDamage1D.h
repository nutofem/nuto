// $Id: Truss.h 627 2013-05-22 07:43:22Z unger3 $
#ifndef BOUNDARYGRADIENTDAMAGE1D_H
#define BOUNDARYGRADIENTDAMAGE1D_H

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/Truss1D.h"

namespace NuTo
{

namespace BoundaryCondition
{

enum eType
{
    NOT_SET,
    NEUMANN_HOMOGENEOUS,            // grad nonlocal eq strain * n = 0
    DIRICHLET_INHOMOGENEOUS,        // nonlocal eq strain = local eq strain
    ROBIN_INHOMOGENEOUS,            // l * grad nonlocal eq strain * n + nonlocal eq strain = local eq strain
    ROBIN_HOMOGENEOUS               // l * grad nonlocal eq strain * n + nonlocal eq strain = 0
};

}  // namespace BoundaryCondition


class StructureBase;
class ConstitutiveTangentLocal1x1;
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
    		Truss* rRealBoundaryElement,
    		int rSurfaceEdge,
    		ElementData::eElementDataType rElementDataType,
    		IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType
    		);

    BoundaryGradientDamage1D(const StructureBase* rStructure,
            Truss* rRealBoundaryElement,
            NodeBase* rSurfaceNode,
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
        return 1;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
        return GetNode(rLocalNodeNumber);
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

    //! @brief calculates the local coordinates of the nodes
    //! @param localCoordinates vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const
    {}

    //! @brief ... interpolate three-dimensional global point coordinates from one-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... one-dimensional local point coordinates
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    void InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const
    {}

    //! @brief calculates the boundary integral of Nt * c * n * B
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rSerivativeShapeFunctions of the ip for all shape functions
    //! @param rNonlocalGradientRadius
    //! @param rFactor multiplication factor (detJ area..)
    //! @param rKkkMod return matrix with detJ * Nt * c * n * B
    void CalculateKkkMod(
            const std::vector<double>& rShapeFunctions,
            const std::vector<double>& rDerivativeShapeFunctions,
            double rNonlocalGradientRadius, double rNormalVector, double rfactor,
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkkMod) const;

    void CalculateKedMod(
            const std::vector<double>& rShapeFunctions,
            const std::vector<double>& rDerivativeShapeFunctions,
            double rNonlocalGradientRadius, double rfactor,
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKedMod) const;


    //! @brief creates constraints according to the boundary condition type
    //! @param rType type of boundary condition, see enum. This is for try out purposes only. If a nice BC type is found, remove this method
    //! @param rStructure nonconst ptr to the structure to create the constraints
    void ApplyConstraints(BoundaryCondition::eType rType, StructureBase* rStructure);

    //! @brief calculates output data fo the elmement
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    NuTo::Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput);

    //! @brief cast the base pointer, otherwise throws an exception
    const NuTo::BoundaryGradientDamage1D* AsBoundaryGradientDamage1D()const override;

    //! @brief cast the base pointer, otherwise throws an exception
    NuTo::BoundaryGradientDamage1D* AsBoundaryGradientDamage1D() override;

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

    // edge number 0.. left, 1.. right
    int mSurfaceEdge;

    BoundaryCondition::eType mBoundaryConditionType;

    // build global row dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int rNumDispDofs, int rNumNonlocalEqStrainDofs) const;

    //! @brief returns the mSurfaceEdge variable for the node rNode
    int CalculateSurfaceEdge(const NodeBase* rNode) const;



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
