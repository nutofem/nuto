#pragma once

#include "nuto/mechanics/elements/ContinuumBoundaryElement.h"
#include <eigen3/Eigen/Dense>

namespace NuTo
{
//! @author Peter Otto, BAM
//! @date September, 2016
//! @brief ... class for contact discretization by mortar method

template<int TDimSlave, int TDimMaster>
class ContinuumContactElement: public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    ContinuumContactElement() = default;
#endif // ENABLE_SERIALIZATION
public:
    ContinuumContactElement(const std::vector<std::pair<const ContinuumElement<TDimSlave> *, int> > &rElementsSlave,
                            Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster> *, int>, Eigen::Dynamic, Eigen::Dynamic> &rElementsMaster,
                            const ConstitutiveBase *rConstitutiveContactLaw,
                            int rContactAlgorithm);

    virtual ~ContinuumContactElement() = default;

    NuTo::eError Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) override;

    NuTo::ConstitutiveOutputMap GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) const;

    void FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian1) const;
    void FillConstitutiveOutputMapHessian2Lumped(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullVector<double> &rHessian2Lumped) const;
    void FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double> &rHessian2) const;

    virtual NuTo::Element::eElementType GetEnumType() const override;

    // abstract functions from element base - not needed for contact element //
    virtual int GetLocalDimension() const override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual NodeBase* GetNode(int rLocalNodeNumber) override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual const NodeBase* GetNode(int rLocalNodeNumber) const override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType) override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual const NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType) const override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual void SetNode(int rLocalNodeNumber, NodeBase* rNode) override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual void ResizeNodes(int rNewNumNodes) override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr) override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual void SetSection(const SectionBase* rSection) {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual const SectionBase* GetSection() const override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}
    virtual const Eigen::VectorXd GetIntegrationPointVolume() const override {throw MechanicsException(__PRETTY_FUNCTION__, "not implemented.");}

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    virtual void CheckElement()
    {
        return;
    }

    void ExtractAllNecessaryDofValues(EvaluateDataContinuumBoundary<TDimSlave>& data, const std::pair<const ContinuumElement<TDimSlave>*, int> &rElementAndSurfaceId);

    void CalculateGlobalRowDofs(BlockFullVector<int> &rGlobalRowDofs) const;

    void CalculateGlobalRowDofsLocal(BlockFullVector<int> &rGlobalRowDofs) const;

    void CalculateGlobalColumnDofs(BlockFullVector<int> &rGlobalDofMapping) const;

    void AddGlobalRowDofsElementMaster(int row, int col) const;

    void GapMatrixMortar(EvaluateDataContinuumBoundary<TDimSlave> &rData,
                         const std::pair<const ContinuumElement<TDimSlave> *, int> &rElementAndSurfaceId);

    void CalculateElementOutputsLocal(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> &rElementOutput,
                                      EvaluateDataContinuumBoundary<TDimSlave>                       &rData,
                                      const std::pair<const ContinuumElement<TDimSlave> *, int>      &rElementAndSurfaceId);

    void CalculateElementOutputsLocalForce(const EvaluateDataContinuumBoundary<TDimSlave> &rData, const std::pair<const ContinuumElement<TDimSlave>*,int> &rElementAndSurfaceId);
    void CalculateElementOutputsLocalForceDerivative(const EvaluateDataContinuumBoundary<TDimSlave> &rData, const std::pair<const ContinuumElement<TDimSlave>*,int> &rElementAndSurfaceId);
    void CalculateElementOutputsLocalGapMatrix(const EvaluateDataContinuumBoundary<TDimSlave> &rData,
                                               const std::pair<const ContinuumElement<TDimSlave>*,int>  &rElementAndSurfaceId,
                                               Eigen::MatrixXd& D, Eigen::MatrixXd& M);

    void CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> &rElementOutput) const;


    void CalculateElementOutputContactForce(BlockFullVector<double>& rInternalGradient) const;

    void CalculateElementOutputContactForceDerivative(BlockFullMatrix<double> &rGapMatrix) const;

    void GetGlobalIntegrationPointCoordinatesAndParameters(int rIpNum,
                                                           Eigen::VectorXd &rCoordinatesIPSlave,
                                                           Eigen::VectorXd &rParamsIPSlave,
                                                           const std::pair<const ContinuumElement<TDimSlave> *, int> &rElementAndSurfaceId) const;

    void ComputeGapMatrix(Eigen::MatrixXd& D, Eigen::MatrixXd& M);

    const ContinuumContactElement<1,1>& AsContinuumContactElement11() const override
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<1,1>.");}

    ContinuumContactElement<1,1>& AsContinuumContactElement11() override
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<1,1>.");}

    const ContinuumContactElement<2,1>& AsContinuumContactElement21() const override
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<1,1>.");}

    ContinuumContactElement<2,1>& AsContinuumContactElement21() override
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumContactElement<1,1>.");}

protected:
    std::vector<std::pair<const ContinuumElement<TDimSlave>*, int> > mElementsSlave;

    //! @brief ... Master elements in the right ordering. For pure IGA structure and contact only 2D elements possible, with the surface id of the contacting surface
    //! in the case of IGA layer ontop of the FEM mesh, there is no need for the surface id => -1
    Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster>*, int>, Eigen::Dynamic, Eigen::Dynamic> mElementsMaster;

    //! @brief ... Matrix containing the knots of all elements, ascending order
    std::vector<Eigen::VectorXd> mKnots;

    std::unordered_map<int, int> mMappingGlobal2LocalDof;

    std::unordered_map<int, int> mMappingGlobal2LocalSlaveNodes;

    int mNumDofs;
    int mNumSlaveDofs;
    int mNumMasterDofs;
    int mNumSlaveNodes;
    bool mDofMappingComputed;
    bool mSlaveNodesMappingComputed;

    const ConstitutiveBase* mConstitutiveContactLaw;

    int mContactType;


    Eigen::VectorXd mContactForce;
    Eigen::VectorXd mMortarGap;

    Eigen::MatrixXd mDerivativeContactForce;

    void FillMappingGlobalLocalDofs();

    void FillMappingGlobalLocalSlaveNodes();
};
} /* namespace NuTo */

