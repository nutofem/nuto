#pragma once

#include "nuto/mechanics/elements/ContinuumBoundaryElement.h"
#include <eigen3/Eigen/Dense>

namespace NuTo
{
//! @author Peter Otto, BAM
//! @date September, 2016
//! @brief ... class for contact discretization by mortar method

template<int TDimSlave, int TDimMaster>
class ContinuumContactElement: public ContinuumBoundaryElement<TDimSlave>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    ContinuumContactElement() = default;
#endif // ENABLE_SERIALIZATION
public:
    ContinuumContactElement(const ContinuumElement<TDimSlave> *rSlaveElement,
                            int rSurfaceId,
                            Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster> *, int>, Eigen::Dynamic, Eigen::Dynamic> &rElementsMaster,
                            double rPenalty, int rContactAlgorithm);

    virtual ~ContinuumContactElement() = default;

    NuTo::eError Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) override;

    NuTo::ConstitutiveOutputMap GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) const override;

    void FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian1) const;
    void FillConstitutiveOutputMapHessian2Lumped(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullVector<double> &rHessian2Lumped) const;
    void FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double> &rHessian2) const;

    NuTo::Element::eElementType GetEnumType() const;

    void CalculateGlobalRowDofs(BlockFullVector<int> &rGlobalRowDofs) const;

    void CalculateGlobalColumnDofs(BlockFullVector<int> &rGlobalDofMapping) const;

    void AddGlobalRowDofsElementMaster(int row, int col) const;

    void GapMatrixMortar(EvaluateDataContinuumBoundary<TDimSlave> &rData, int rTheIP);

    void CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> &rElementOutput,
                                 EvaluateDataContinuumBoundary<TDimSlave>                       &rData,
                                 const ConstitutiveInputMap                                     &constitutiveInput,
                                 const ConstitutiveOutputMap                                    &constitutiveOutput) const;


    void CalculateElementOutputContactForce(BlockFullVector<double>& rInternalGradient,
                                            EvaluateDataContinuumBoundary<TDimSlave> &rData,
                                            const ConstitutiveInputMap& constitutiveInput,
                                            const ConstitutiveOutputMap& constitutiveOutput) const;

    void CalculateElementOutputContactForceDerivative(BlockFullMatrix<double> &rGapMatrix,
                                                      EvaluateDataContinuumBoundary<TDimSlave> &rData,
                                                      const ConstitutiveInputMap &constitutiveInput,
                                                      const ConstitutiveOutputMap &constitutiveOutput) const;

    void GetGlobalIntegrationPointCoordinatesAndParameters(int rIpNum, Eigen::VectorXd &rCoordinatesIPSlave, Eigen::VectorXd &rParamsIPSlave) const;

    const Eigen::Vector3d GetGlobalIntegrationPointCoordinates(int rIpNum) const override;


protected:

    //! @brief ... Master elements in the right ordering. For pure IGA structure and contact only 2D elements possible, with the surface id of the contacting surface
    //! in the case of IGA layer ontop of the FEM mesh, there is no need for the surface id => -1
    Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster>*, int>, Eigen::Dynamic, Eigen::Dynamic> mElementsMaster;

    //! @brief ... Matrix containing the knots of all elements, ascending order
    std::vector<Eigen::VectorXd> mKnots;

    std::unordered_map<int, int> mMappingGlobal2LocalDof;

    int mNumDofs;

    double penaltyP;

    int mContactType;

    void FillMappingGlobalLocal();
};
} /* namespace NuTo */

