#pragma once

#include "ContinuumBoundaryElement.h"
#include <eigen3/Eigen/Dense>


namespace NuTo
{
//! @author Peter Otto, BAM
//! @date September, 2016
//! @brief ... class for contact discretization by mortar method

template <int TDimSlave, int TDimMaster>
class ContinuumContactElement : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;

protected:
    ContinuumContactElement() = default;
#endif // ENABLE_SERIALIZATION
public:
    ContinuumContactElement(const std::vector<std::pair<const ContinuumElement<TDimSlave>*, int>>& rElementsSlave,
                            Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster>*, int>, Eigen::Dynamic,
                                          Eigen::Dynamic>& rElementsMaster,
                            const ConstitutiveBase* rConstitutiveContactLaw, int rContactAlgorithm,
                            const std::function<bool(int, int)>& IsNodeOnSurface);

    virtual ~ContinuumContactElement()
    {
    }

    NuTo::eError Evaluate(const ConstitutiveInputMap& rInput,
                          std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) override;

    NuTo::ConstitutiveOutputMap
    GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) const;

    void FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap& rConstitutiveOutput,
                                           BlockFullMatrix<double>& rHessian1) const;
    void FillConstitutiveOutputMapHessian2Lumped(ConstitutiveOutputMap& rConstitutiveOutput,
                                                 BlockFullVector<double>& rHessian2Lumped) const;
    void FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap& rConstitutiveOutput,
                                           BlockFullMatrix<double>& rHessian2) const;

    virtual NuTo::Element::eElementType GetEnumType() const override;

    // abstract functions from element base - not needed for contact element //
    virtual int GetLocalDimension() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual NodeBase* GetNode(int rLocalNodeNumber) override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual const NodeBase* GetNode(int rLocalNodeNumber) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType) override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual const NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual void SetNode(int rLocalNodeNumber, NodeBase* rNode) override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual void ResizeNodes(int rNewNumNodes) override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr) override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual void SetSection(const SectionBase* rSection)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual const SectionBase* GetSection() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }
    virtual const Eigen::VectorXd GetIntegrationPointVolume() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");
    }

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element
    //! length/area/volum is negative)
    virtual void CheckElement()
    {
        return;
    }

    void ExtractAllNecessaryDofValues(EvaluateDataContinuumBoundary<TDimSlave>& data,
                                      const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId);

    void CalculateGlobalRowDofs(BlockFullVector<int>& rGlobalRowDofs) const;

    void CalculateGlobalRowDofsLocal(BlockFullVector<int>& rGlobalRowDofs) const;

    void CalculateGlobalColumnDofs(BlockFullVector<int>& rGlobalDofMapping) const;

    void AddGlobalRowDofsElementMaster(int row, int col) const;

    void GapMatrixMortarContact(EvaluateDataContinuumBoundary<TDimSlave>& rData,
                                const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId,
                                bool rAssembleForce, bool rAssembleDerivative);

    void GapMatrixMortarTying(const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId,
                              Eigen::MatrixXd& D, Eigen::MatrixXd& M);

    const ContinuumElementIGA<TDimMaster>*
    Projection(const Eigen::VectorXd& coordinatesIPSlave,
               const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId,
               Eigen::VectorXd& rProjectionVector, Eigen::VectorXd& rParameterMinMaster);

    void AssembleMortarTypeContact(const EvaluateDataContinuumBoundary<TDimSlave>& rData,
                                   const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId,
                                   bool rAssembleGapVector = true);

    void AssembleNonMortarTypeContact(const EvaluateDataContinuumBoundary<TDimSlave>& rData,
                                      const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId,
                                      bool rAssembleForce, bool rAssembleDerivative);

    void CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput);

    void CalculateElementOutputContactForce(BlockFullVector<double>& rInternalGradient);

    void CalculateElementOutputContactForceDerivative(BlockFullMatrix<double>& rGapMatrix);

    void CalculateElementOutputContactForceDerivative(Eigen::MatrixXd& der);

    void GetGlobalIntegrationPointCoordinatesAndParameters(
            int rIpNum, Eigen::VectorXd& rCoordinatesIPSlave, Eigen::VectorXd& rParamsIPSlave,
            const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId) const;

    void ComputeIndicesForElementAssemblyMeshTying(
            const std::pair<const ContinuumElement<TDimSlave>*, int>& rSlaveElementAndSurfaceId,
            const ContinuumElement<TDimMaster>* rMasterElement, Eigen::VectorXi& indicesNodesSlave,
            Eigen::VectorXi& indicesNodesMaster);

    void ComputeMeshTyingMatrix(Eigen::MatrixXd& D, Eigen::MatrixXd& M,
                                std::unordered_map<int, int>& mappingGlobal2LocalSlaveNode,
                                std::unordered_map<int, int>& mappingGlobal2LocalMasterNode);

    Eigen::VectorXd GetContactPressure(std::unordered_map<int, int>& rMappingGlobal2LocalSlaveNodes);

    double CalculateContactForce();

    std::unordered_map<int, int> GetDofMapping()
    {
        return mMappingGlobal2LocalDof;
    }

    const ContinuumContactElement<1, 1>& AsContinuumContactElement11() const override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<1,1>.");
    }

    ContinuumContactElement<1, 1>& AsContinuumContactElement11() override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<1,1>.");
    }

    const ContinuumContactElement<2, 1>& AsContinuumContactElement21() const override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<2,1>.");
    }

    ContinuumContactElement<2, 1>& AsContinuumContactElement21() override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<2,1>.");
    }

    const ContinuumContactElement<2, 2>& AsContinuumContactElement22() const override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<2,2>.");
    }

    ContinuumContactElement<2, 2>& AsContinuumContactElement22() override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<2,2>.");
    }

    ContinuumContactElement<3, 2>& AsContinuumContactElement32() override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<3,2>.");
    }


protected:
    std::vector<std::pair<const ContinuumElement<TDimSlave>*, int>> mElementsSlave;

    //! @brief ... Master elements in the right ordering. For pure IGA structure and contact only 2D elements possible,
    //! with the surface id of the contacting surface
    //! in the case of IGA layer ontop of the FEM mesh, there is no need for the surface id => -1
    Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster>*, int>, Eigen::Dynamic, Eigen::Dynamic>
            mElementsMaster;

    //! @brief ... Matrix containing the knots of all elements, ascending order
    std::vector<Eigen::VectorXd> mKnots;

    std::unordered_map<int, int> mMappingGlobal2LocalDof;
    std::unordered_map<int, int> mMappingGlobal2LocalSlaveNodes;
    std::unordered_map<int, int> mMappingGlobal2LocalMasterNodes;

    int mNumDofs;
    int mNumSlaveDofs;
    int mNumMasterDofs;

    int mNumSlaveNodes;
    int mNumMasterNodes;

    bool mDofMappingComputed;
    bool mSlaveNodesMappingComputed;
    bool mMasterNodesMappingComputed;


    const ConstitutiveBase* mConstitutiveContactLaw;

    int mContactType;

    // *** the result of the contact algorithm *** //

    Eigen::MatrixXd mGapMatrixMeshTying; // both mortar and non-mortar
    bool mGapMatrixMeshTyingAssembled;

    Eigen::MatrixXd mGapMatrix; // both mortar and non-mortar
    bool mGapMatrixAssembled;

    Eigen::VectorXd mMortarGlobalGapVector; // contact type = 0: mortar
    bool mMortarGlobalGapVectorAssembled;

    Eigen::VectorXd mSlaveShapeFunctionsWeight; // int_{\Gamma_c}R_A^s\;d\Gamma
    bool mSlaveShapeFunctionsWeightAssembled;

    Eigen::VectorXd mGlobalNodalPressure; // contact type = 1: nonmortar
    bool mGlobalNodalPressureAssembled;

    Eigen::MatrixXd mGapMatrixPenalty; // contact type = 1: nonmortar
    bool mGapMatrixPenaltyAssembled;

    //    Eigen::MatrixXd mGapMatrixTemp;
    //    Eigen::VectorXd mGapVectorTemp;
    //    Eigen::VectorXd mPressureVectorTemp;
    //    Eigen::VectorXd mWeightsVectorTemp;

    std::function<bool(int, int)> mIsNodeOnSurface;

    // *** the assembly takes place here, so mapping needed local - global dof numbers *** //

    void FillMappingGlobalLocalDofs();
    void FillMappingGlobalLocalSlaveNodes();
    void FillMappingGlobalLocalMasterNodes();
};
} /* namespace NuTo */
