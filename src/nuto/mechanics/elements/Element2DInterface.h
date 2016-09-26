//============================================================================
// Name        : Element2DInterface.cpp
// Author      : Philip Huschke
// Version     : 26 Aug 2015
// Copyright   :
// Description : Element formulation for the interface element proposed by Goodman et al.
//============================================================================

#pragma once

#include "nuto/mechanics/elements/ElementBase.h"

namespace NuTo
{
class ElementOutputIpData;
template <typename T> class BlockFullVector;
template <typename T> class BlockFullMatrix;

struct EvaluateData
{
    std::map<Node::eDof, Eigen::VectorXd> mNodalValues;
    std::map<Node::eDof, Eigen::MatrixXd> mMatrixB;
    std::map<Node::eDof, const Eigen::MatrixXd*> mMatrixN;

    double mDetJacobian = 0;
    double mDetJxWeightIPxSection = 0;
};


class Element2DInterface: public ElementBase
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    Element2DInterface(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType);


    Element2DInterface(const Element2DInterface& ) = default;
    Element2DInterface(      Element2DInterface&&) = default;
    virtual ~Element2DInterface() = default;




    //! @brief calculates output data for the element
    //! @param rInput ... constitutive input map for the constitutive law
    //! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    eError Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput) override;


    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override;

    //! @brief returns the local dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! @return local dimension
    int GetLocalDimension() const override
    {
        return 1;
    }


    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    Constitutive::StaticData::Component* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber) override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber) const override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType) override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType) const override;

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode) override;

    //! @brief resizes the node vector
    //! @param rNewNumNodes new number of nodes
    void ResizeNodes(int rNewNumNodes) override;

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr) override;

    virtual Eigen::VectorXd ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const override;

    //! @brief sets the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @param rSection pointer to section
    void SetSection(const SectionBase* rSection) override
    {
        mSection = rSection;
    }

    //! @brief returns a pointer tclass ElementOutputIpData;o the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @return pointer to section
    const SectionBase* GetSection() const override
    {
        return mSection;
    }

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @return rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    const Eigen::VectorXd GetIntegrationPointVolume() const override;




protected:

    std::vector<NodeBase*> mNodes;
    const SectionBase *mSection;
    Eigen::MatrixXd mTransformationMatrix;


    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override;

private:
    ConstitutiveOutputMap GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase> >& rElementOutput);
    ConstitutiveInputMap GetConstitutiveInputMap(const ConstitutiveOutputMap& rConstitutiveOutput) const;
    void CalculateConstitutiveInputs(const ConstitutiveInputMap& rConstitutiveInput, EvaluateData& rData);

    virtual void FillConstitutiveOutputMapInternalGradient(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullVector<double>& rInternalGradient) const;
    virtual void FillConstitutiveOutputMapHessian0(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian0) const;
    virtual void FillConstitutiveOutputMapIpData(ConstitutiveOutputMap& rConstitutiveOutput, ElementOutputIpData& rIpData) const;

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    void CalculateGlobalRowDofs(BlockFullVector<int>& rGlobalRowDofs) const;

    void CalculateElementOutputs(
            std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
            EvaluateData& rData, int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const;

    virtual void CalculateElementOutputInternalGradient(BlockFullVector<double>& rInternalGradient,
            EvaluateData& rData, int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const;
    virtual void CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0,
            EvaluateData& rData, int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const;
    virtual void CalculateElementOutputIpData(              ElementOutputIpData&     rIpData,           EvaluateData& rData, int rTheIP) const;

    //! @brief calculates the rotation matirx based on the orientation of the element
    Eigen::MatrixXd CalculateRotationMatrix();

    //! @brief calculates the transformation matrix, which contains the rotation matrices on the diagonal
    Eigen::MatrixXd CalculateTransformationMatrix(unsigned int rGlobalDimension, unsigned int rNumberOfNodes);

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
            std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const override;

    void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) override;
#endif // ENABLE_VISUALIZE
};

} /* namespace NuTo */





















