/*
 * ContinuumBoundaryElementBase.cpp
 *
 *  Created on: 5 Mar 2016
 *      Author: vhirtham
 */
#include <iostream>

#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/elements/ContinuumBoundaryElement.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputBase.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/EvaluateDataContinuumBoundary.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/sections/Section.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

using namespace NuTo;

template <int TDim>
NuTo::ContinuumBoundaryElement<TDim>::ContinuumBoundaryElement(const ContinuumElement<TDim>& rBaseElement,
                                                               const IntegrationTypeBase& integrationType,
                                                               int rSurfaceId)
    : ElementBase::ElementBase(rBaseElement.GetInterpolationType(), integrationType)
    , mBaseElement(rBaseElement)
    , mSurfaceId(rSurfaceId)
    , mAlphaUserDefined(-1)
{
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::Evaluate(
        const ConstitutiveInputMap& rInput,
        std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput)
{
    EvaluateDataContinuumBoundary<TDim> data;
    ExtractAllNecessaryDofValues(data);

    auto constitutiveOutput = GetConstitutiveOutputMap(rElementOutput);
    auto constitutiveInput = GetConstitutiveInputMap(constitutiveOutput);
    constitutiveInput.Merge(rInput);

    if (TDim == 2)
        AddPlaneStateToInput(constitutiveInput);

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        CalculateNMatrixBMatrixDetJacobian(data, theIP);
        CalculateConstitutiveInputs(constitutiveInput, data);

        EvaluateConstitutiveLaw<TDim>(constitutiveInput, constitutiveOutput, theIP);
        CalculateElementOutputs(rElementOutput, data, theIP, constitutiveInput, constitutiveOutput);
    }
}


template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::ExtractAllNecessaryDofValues(EvaluateDataContinuumBoundary<TDim>& rData)
{
    // needs optimization,
    // not all dofs might be needed...

    const std::set<Node::eDof>& dofs = mInterpolationType->GetDofs();
    for (auto dof : dofs)
        if (mInterpolationType->IsConstitutiveInput(dof))
            rData.mNodalValues[dof] = mBaseElement.ExtractNodeValues(0, dof);

    rData.mNodalValues[Node::eDof::COORDINATES] = mBaseElement.ExtractNodeValues(0, Node::eDof::COORDINATES);

    for (auto dof : dofs)
        if (mInterpolationType->IsConstitutiveInput(dof))
            if (GetNode(0)->GetNumTimeDerivatives(dof) >= 1)
                rData.mNodalValues_dt1[dof] = mBaseElement.ExtractNodeValues(1, dof);
}


template <int TDim>
NuTo::ConstitutiveOutputMap NuTo::ContinuumBoundaryElement<TDim>::GetConstitutiveOutputMap(
        std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) const
{
    ConstitutiveOutputMap constitutiveOutput;

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case ElementEnum::eOutput::INTERNAL_GRADIENT:
            FillConstitutiveOutputMapInternalGradient(constitutiveOutput, it.second->GetBlockFullVectorDouble());
            break;

        case ElementEnum::eOutput::HESSIAN_0_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian0(constitutiveOutput, it.second->GetBlockFullMatrixDouble());
            break;

        case ElementEnum::eOutput::HESSIAN_1_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian1(constitutiveOutput, it.second->GetBlockFullMatrixDouble());
            break;

        case ElementEnum::eOutput::HESSIAN_2_TIME_DERIVATIVE:
            throw Exception(__PRETTY_FUNCTION__,
                                     "Case not handled! IMPORTANT: Everything must be set to zero here!!!");
            break;

        case ElementEnum::eOutput::UPDATE_STATIC_DATA:
            constitutiveOutput[Constitutive::eOutput::UPDATE_STATIC_DATA] = nullptr;
            break;

        case ElementEnum::eOutput::UPDATE_TMP_STATIC_DATA:
            constitutiveOutput[Constitutive::eOutput::UPDATE_TMP_STATIC_DATA] = nullptr;
            break;

        case ElementEnum::eOutput::IP_DATA:
            FillConstitutiveOutputMapIpData(constitutiveOutput, it.second->GetIpData());
            break;

        case ElementEnum::eOutput::GLOBAL_ROW_DOF:
            mBaseElement.CalculateGlobalRowDofs(it.second->GetBlockFullVectorInt());
            break;

        case ElementEnum::eOutput::GLOBAL_COLUMN_DOF:
            mBaseElement.CalculateGlobalColumnDofs(it.second->GetBlockFullVectorInt());
            break;

        default:
            throw Exception(__PRETTY_FUNCTION__, "element  output not implemented.");
        }
    }

    // allocate the objects for the output data
    for (auto& outputs : constitutiveOutput)
    {
        outputs.second = ConstitutiveIOBase::makeConstitutiveIO<TDim>(outputs.first);
    }
    return constitutiveOutput;
}


template <int TDim>
NuTo::ConstitutiveInputMap
NuTo::ContinuumBoundaryElement<TDim>::GetConstitutiveInputMap(const ConstitutiveOutputMap& rConstitutiveOutput) const
{
    ConstitutiveInputMap constitutiveInput = GetConstitutiveLaw(0).GetConstitutiveInputs(rConstitutiveOutput);

    for (auto& itInput : constitutiveInput)
    {
        itInput.second = ConstitutiveIOBase::makeConstitutiveIO<TDim>(itInput.first);
    }
    return constitutiveInput;
}


template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateNMatrixBMatrixDetJacobian(
        EvaluateDataContinuumBoundary<TDim>& rData, int rTheIP) const
{

    const InterpolationBase& interpolationTypeCoords = mInterpolationType->Get(Node::eDof::COORDINATES);

    Eigen::Matrix<double, TDim - 1, 1> ipCoordsSurface = CalculateIPCoordinatesSurface(rTheIP);
    Eigen::Matrix<double, TDim, 1> ipCoordsNatural =
            interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, mSurfaceId);

    // #######################################
    // ##  Calculate the surface jacobian
    // ## = || [dX / dXi] * [dXi / dAlpha] ||
    // #######################################
    Eigen::MatrixXd derivativeShapeFunctionsNatural =
            interpolationTypeCoords.DerivativeShapeFunctionsNatural(ipCoordsNatural);
    const Eigen::Matrix<double, TDim, TDim> jacobian = mBaseElement.CalculateJacobian(
            derivativeShapeFunctionsNatural, rData.mNodalValues[Node::eDof::COORDINATES]); // = [dX / dXi]

    const Eigen::MatrixXd derivativeNaturalSurfaceCoordinates =
            interpolationTypeCoords.CalculateDerivativeNaturalSurfaceCoordinates(ipCoordsSurface,
                                                                                 mSurfaceId); // = [dXi / dAlpha]
    rData.mDetJacobian = (jacobian * derivativeNaturalSurfaceCoordinates).norm();

    if (rData.mDetJacobian == 0)
    {
        throw Exception(__PRETTY_FUNCTION__, "Determinant of the Jacobian is zero, no inversion possible.");
    }

    const Eigen::Matrix<double, TDim, TDim> invJacobian = jacobian.inverse();


    for (auto dof : mInterpolationType->GetDofs())
    {
        if (dof == Node::eDof::COORDINATES)
            continue;
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        rData.mN[dof] = interpolationType.MatrixN(ipCoordsNatural);

        rData.mB[dof] = mBaseElement.CalculateMatrixB(
                dof, interpolationType.DerivativeShapeFunctionsNatural(ipCoordsNatural), invJacobian);
    }
}


template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateConstitutiveInputs(const ConstitutiveInputMap& rConstitutiveInput,
                                                                       EvaluateDataContinuumBoundary<TDim>& rData)
{
    constexpr int voigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    for (auto& it : rConstitutiveInput)
    {
        switch (it.first)
        {
        case Constitutive::eInput::ENGINEERING_STRAIN:
        {
            auto& strain = *static_cast<ConstitutiveVector<voigtDim>*>(it.second.get());
            strain.AsVector() =
                    rData.mB.at(Node::eDof::DISPLACEMENTS) * rData.mNodalValues.at(Node::eDof::DISPLACEMENTS);
            break;
        }
        case Constitutive::eInput::NONLOCAL_EQ_STRAIN:
        {
            auto& nonLocalEqStrain = *static_cast<ConstitutiveScalar*>(it.second.get());
            nonLocalEqStrain.AsScalar() =
                    rData.mN.at(Node::eDof::NONLOCALEQSTRAIN) * rData.mNodalValues.at(Node::eDof::NONLOCALEQSTRAIN);
            break;
        }
        case Constitutive::eInput::RELATIVE_HUMIDITY:
        {
            auto& relativeHumidity = *static_cast<ConstitutiveScalar*>(it.second.get());
            relativeHumidity.AsScalar() =
                    rData.mN.at(Node::eDof::RELATIVEHUMIDITY) * rData.mNodalValues.at(Node::eDof::RELATIVEHUMIDITY);
            break;
        }
        case Constitutive::eInput::RELATIVE_HUMIDITY_BOUNDARY:
        {
            auto& relativeHumidityBoundary = *static_cast<ConstitutiveScalar*>(it.second.get());
            relativeHumidityBoundary.AsScalar() = this->GetBoundaryControlNode()->Get(Node::eDof::RELATIVEHUMIDITY);
            break;
        }
        case Constitutive::eInput::WATER_VOLUME_FRACTION:
        {
            auto& waterVolumeFraction = *static_cast<ConstitutiveScalar*>(it.second.get());
            waterVolumeFraction.AsScalar() = rData.mN.at(Node::eDof::WATERVOLUMEFRACTION) *
                                             rData.mNodalValues.at(Node::eDof::WATERVOLUMEFRACTION);
            break;
        }
        case Constitutive::eInput::TIME:
        case Constitutive::eInput::TIME_STEP:
        case Constitutive::eInput::CALCULATE_STATIC_DATA:
        case Constitutive::eInput::PLANE_STATE:
            break;

        default:
            throw Exception(__PRETTY_FUNCTION__, "Constitutive input for " +
                                                                  Constitutive::InputToString(it.first) +
                                                                  " not implemented.");
        }
    }
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputs(
        std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
        EvaluateDataContinuumBoundary<TDim>& rData, int rTheIP, const ConstitutiveInputMap& constitutiveInput,
        const ConstitutiveOutputMap& constitutiveOutput) const
{
    rData.mDetJxWeightIPxSection =
            CalculateDetJxWeightIPxSection(rData.mDetJacobian, rTheIP); // formerly known as "factor"

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case ElementEnum::eOutput::INTERNAL_GRADIENT:
            UpdateAlphaGradientDamage(rData, constitutiveInput, constitutiveOutput);
            CalculateElementOutputInternalGradient(it.second->GetBlockFullVectorDouble(), rData, constitutiveInput,
                                                   constitutiveOutput, rTheIP);
            break;

        case ElementEnum::eOutput::HESSIAN_0_TIME_DERIVATIVE:
            UpdateAlphaGradientDamage(rData, constitutiveInput, constitutiveOutput);
            CalculateElementOutputHessian0(it.second->GetBlockFullMatrixDouble(), rData, constitutiveOutput, rTheIP);
            break;

        case ElementEnum::eOutput::HESSIAN_1_TIME_DERIVATIVE:
            break;

        case ElementEnum::eOutput::HESSIAN_2_TIME_DERIVATIVE:
            break;

        case ElementEnum::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
            break;

        case ElementEnum::eOutput::UPDATE_STATIC_DATA:
        case ElementEnum::eOutput::UPDATE_TMP_STATIC_DATA:
            break;
        case ElementEnum::eOutput::IP_DATA:
            CalculateElementOutputIpData(it.second->GetIpData(), constitutiveOutput, rTheIP);
            break;
        case ElementEnum::eOutput::GLOBAL_ROW_DOF:
        case ElementEnum::eOutput::GLOBAL_COLUMN_DOF:
            break;
        default:
            throw Exception(__PRETTY_FUNCTION__, "element output not implemented.");
        }
    }
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::UpdateAlphaGradientDamage(
        NuTo::EvaluateDataContinuumBoundary<TDim>& rData, const NuTo::ConstitutiveInputMap&,
        const NuTo::ConstitutiveOutputMap&) const
{
    if (this->mInterpolationType->GetActiveDofs().find(NuTo::Node::eDof::NONLOCALEQSTRAIN) !=
        this->mInterpolationType->GetActiveDofs().end())
    {
        if (this->mAlphaUserDefined == -1) // just to put it somewhere...
            rData.m1DivAlpha = 1. / this->CalculateAlpha();
        else
            rData.m1DivAlpha = 1. / this->mAlphaUserDefined;

        // const auto &localEqStrain =
        //    *static_cast<NuTo::ConstitutiveScalar
        //    *>(rConstitutiveOutput.at(NuTo::Constitutive::eOutput::LOCAL_EQ_STRAIN).get());
        // const auto &nonlocalEqStrain =
        //    *static_cast<NuTo::ConstitutiveScalar
        //    *>(rConstitutiveInput.at(NuTo::Constitutive::eInput::NONLOCAL_EQ_STRAIN).get());
        //
        // double gradEeqN = localEqStrain[0] - nonlocalEqStrain[0];
        // if (gradEeqN > 0)
        //    rData.m1DivAlpha = 0;
    }
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputInternalGradient(
        BlockFullVector<double>& rInternalGradient, EvaluateDataContinuumBoundary<TDim>& rData,
        const ConstitutiveInputMap& constitutiveInput, const ConstitutiveOutputMap& constitutiveOutput,
        int) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
            break;

        case Node::eDof::NONLOCALEQSTRAIN:
        {
            const auto& localEqStrain = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::LOCAL_EQ_STRAIN).get());
            const auto& nonlocalEqStrain = *static_cast<ConstitutiveScalar*>(
                    constitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN).get());

            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * rData.m1DivAlpha *
                                         rData.mN.at(Node::eDof::NONLOCALEQSTRAIN).transpose() *
                                         (nonlocalEqStrain[0] - localEqStrain[0]);
            break;
        }

        case Node::eDof::RELATIVEHUMIDITY:
        {
            const auto& internalGradientRH_Boundary_N = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N).get());
            rInternalGradient[dofRow] +=
                    rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * internalGradientRH_Boundary_N;
            break;
        }
        case Node::eDof::WATERVOLUMEFRACTION:
        {
            const auto& internalGradientWV_Boundary_N = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N)
                            .get());
            rInternalGradient[dofRow] +=
                    rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * internalGradientWV_Boundary_N;
            break;
        }
        default:
            throw Exception(__PRETTY_FUNCTION__, "Element output INTERNAL_GRADIENT for " +
                                                                  Node::DofToString(dofRow) + " not implemented.");
        }
    }
}


template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputHessian0(
        BlockFullMatrix<double>& rHessian0, EvaluateDataContinuumBoundary<TDim>& rData,
        const ConstitutiveOutputMap& constitutiveOutput, int rTheIP) const
{
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            if (!GetConstitutiveLaw(rTheIP).CheckDofCombinationComputable(dofRow, dofCol, 0))
                continue;
            auto& hessian0 = rHessian0(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN):
                break;

            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::DISPLACEMENTS):
            {
                const auto& tangentLocalEqStrainStrain = *static_cast<ConstitutiveVector<VoigtDim>*>(
                        constitutiveOutput.at(Constitutive::eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN).get());
                hessian0 -= rData.mDetJxWeightIPxSection * rData.m1DivAlpha * rData.mN.at(dofRow).transpose() *
                            tangentLocalEqStrainStrain.transpose() * rData.mB.at(dofCol);
                break;
            }
            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::NONLOCALEQSTRAIN):
            {
                hessian0 += rData.mN.at(dofRow).transpose() * rData.mN.at(dofRow) * rData.mDetJxWeightIPxSection *
                            rData.m1DivAlpha;
                break;
            }
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
            {
                const auto& internalGradientRH_dRH_Boundary_NN_H0 = *static_cast<ConstitutiveScalar*>(
                        constitutiveOutput.at(Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0).get());
                hessian0 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() *
                            internalGradientRH_dRH_Boundary_NN_H0 * rData.mN.at(dofCol);
                break;
            }
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
            {
                const auto& internalGradientWV_dWV_Boundary_NN_H0 = *static_cast<ConstitutiveScalar*>(
                        constitutiveOutput.at(Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0).get());
                hessian0 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() *
                            internalGradientWV_dWV_Boundary_NN_H0 * rData.mN.at(dofCol);
                break;
            }
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
                break;
            default:
            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::DISPLACEMENTS):
                continue;
                throw Exception(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for "
                                                              "(" + Node::DofToString(dofRow) +
                                                                      "," + Node::DofToString(dofCol) +
                                                                      ") not implemented.");
            }
        }
    }
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputIpData(ElementOutputIpData& rIpData,
                                                                        const ConstitutiveOutputMap& constitutiveOutput,
                                                                        int rTheIP) const
{
    for (auto& it :
         rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    { // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRAIN:
            it.second.col(rTheIP) = *static_cast<EngineeringStrain<TDim>*>(
                    constitutiveOutput.at(Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE).get());
            break;
        case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS:
            it.second.col(rTheIP) = *static_cast<EngineeringStress<TDim>*>(
                    constitutiveOutput.at(Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE).get());
            break;
        case NuTo::IpData::eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN:
            it.second.col(rTheIP) = *static_cast<EngineeringStrain<TDim>*>(
                    constitutiveOutput.at(Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE).get());
            break;
        case NuTo::IpData::eIpStaticDataType::DAMAGE:
            it.second.col(rTheIP) =
                    *static_cast<ConstitutiveScalar*>(constitutiveOutput.at(Constitutive::eOutput::DAMAGE).get());
            break;
        case NuTo::IpData::eIpStaticDataType::EXTRAPOLATION_ERROR:
            it.second.col(rTheIP) = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::EXTRAPOLATION_ERROR).get());
            break;
        case NuTo::IpData::eIpStaticDataType::LOCAL_EQ_STRAIN:
            it.second.col(rTheIP) = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::LOCAL_EQ_STRAIN).get());
            break;
        default:
            throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "] Ip data not implemented.");
        }
    }
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapInternalGradient(
        ConstitutiveOutputMap& rConstitutiveOutput, BlockFullVector<double>& rInternalGradient) const
{
    using namespace NuTo::Constitutive;
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        rInternalGradient[dofRow].resize(mInterpolationType->Get(dofRow).GetNumDofs());
        rInternalGradient[dofRow].setZero();
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
            break;

        case Node::eDof::NONLOCALEQSTRAIN:
            rConstitutiveOutput[eOutput::LOCAL_EQ_STRAIN] =
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::LOCAL_EQ_STRAIN);
            break;

        case Node::eDof::RELATIVEHUMIDITY:
            rConstitutiveOutput[eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N] =
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(
                            eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N);
            break;

        case Node::eDof::WATERVOLUMEFRACTION:
            rConstitutiveOutput[eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N] =
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(
                            eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N);
            break;

        default:
            throw Exception(__PRETTY_FUNCTION__, "Constitutive output INTERNAL_GRADIENT for " +
                                                                  Node::DofToString(dofRow) + " not implemented.");
        }
    }
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapHessian0(ConstitutiveOutputMap& rConstitutiveOutput,
                                                                             BlockFullMatrix<double>& rHessian0) const
{
    using namespace NuTo::Constitutive;
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            Eigen::MatrixXd& dofSubMatrix = rHessian0(dofRow, dofCol);
            dofSubMatrix.resize(mInterpolationType->Get(dofRow).GetNumDofs(),
                                mInterpolationType->Get(dofCol).GetNumDofs());
            dofSubMatrix.setZero();
            if (!GetConstitutiveLaw(0).CheckDofCombinationComputable(dofRow, dofCol, 0))
                continue;

            switch (Node::CombineDofs(dofRow, dofCol))
            {

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN):
                break;

            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::DISPLACEMENTS):
                rConstitutiveOutput[eOutput::LOCAL_EQ_STRAIN] =
                        ConstitutiveIOBase::makeConstitutiveIO<TDim>(Constitutive::eOutput::LOCAL_EQ_STRAIN);
                rConstitutiveOutput[eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN] =
                        ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN);
                break;

            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::NONLOCALEQSTRAIN):
                rConstitutiveOutput[eOutput::LOCAL_EQ_STRAIN] =
                        ConstitutiveIOBase::makeConstitutiveIO<TDim>(Constitutive::eOutput::LOCAL_EQ_STRAIN);
                break;

            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
                rConstitutiveOutput[eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0] =
                        ConstitutiveIOBase::makeConstitutiveIO<TDim>(
                                eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0);
                break;

            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
                rConstitutiveOutput[eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0] =
                        ConstitutiveIOBase::makeConstitutiveIO<TDim>(
                                eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0);
                break;


            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
                continue;
            default:
                throw Exception(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_0_TIME_DERIVATIVE for "
                                                              "(" + Node::DofToString(dofRow) +
                                                                      "," + Node::DofToString(dofCol) +
                                                                      ") not implemented.");
            }
        }
    }
}


template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap&,
                                                                             BlockFullMatrix<double>& rHessian0) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            Eigen::MatrixXd& dofSubMatrix = rHessian0(dofRow, dofCol);
            dofSubMatrix.resize(mInterpolationType->Get(dofRow).GetNumDofs(),
                                mInterpolationType->Get(dofCol).GetNumDofs());
            dofSubMatrix.setZero();
            if (!GetConstitutiveLaw(0).CheckDofCombinationComputable(dofRow, dofCol, 1))
                continue;

            switch (Node::CombineDofs(dofRow, dofCol))
            {


            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
            //            case Node::CombineDofs(Node::eDof::DISPLACEMENTS,         Node::eDof::DISPLACEMENTS):
            //            case Node::CombineDofs(Node::eDof::DISPLACEMENTS,         Node::eDof::RELATIVEHUMIDITY):
            //            case Node::CombineDofs(Node::eDof::DISPLACEMENTS,         Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
            //            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY,      Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
                //            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION,   Node::eDof::DISPLACEMENTS):
                //            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION,   Node::eDof::RELATIVEHUMIDITY):
                continue;
            default:
                throw Exception(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_1_TIME_DERIVATIVE for "
                                                              "(" + Node::DofToString(dofRow) +
                                                                      "," + Node::DofToString(dofCol) +
                                                                      ") not implemented.");
            }
        }
    }
}


// template<int TDim>
// void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap&
// rConstitutiveOutput, BlockFullMatrix<double>& rHessian2, EvaluateDataContinuum<TDim> &rData) const
//{
//    for (auto dofRow : mInterpolationType->GetActiveDofs())
//    {
//        for (auto dofCol : mInterpolationType->GetActiveDofs())
//        {
//            Eigen::MatrixXd& dofSubMatrix = rHessian2(dofRow, dofCol);
//            dofSubMatrix.Resize(mInterpolationType->Get(dofRow).GetNumDofs(),
//            mInterpolationType->Get(dofCol).GetNumDofs());
//            dofSubMatrix.setZero();
//            if(!GetConstitutiveLaw(0)->CheckDofCombinationComputable(dofRow,dofCol,2))
//                continue;
//            switch (Node::CombineDofs(dofRow, dofCol))
//            {
//            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
//                break;
//            default:
//                throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "] Constitutive output
//                HESSIAN_2_TIME_DERIVATIVE for "
//                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
//            }
//        }
//    }
//}
//
template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapIpData(ConstitutiveOutputMap& rConstitutiveOutput,
                                                                           ElementOutputIpData& rIpData) const
{
    using namespace NuTo::Constitutive;
    for (auto& it :
         rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    { // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRAIN:
            it.second.resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[eOutput::ENGINEERING_STRAIN_VISUALIZE] =
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::ENGINEERING_STRAIN_VISUALIZE);
            break;
        case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS:
            it.second.resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[eOutput::ENGINEERING_STRESS_VISUALIZE] =
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::ENGINEERING_STRESS_VISUALIZE);
            break;
        case NuTo::IpData::eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN:
            it.second.resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE] =
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE);
            break;
        case NuTo::IpData::eIpStaticDataType::DAMAGE:
            it.second.resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[eOutput::DAMAGE] = ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::DAMAGE);
            break;
        case NuTo::IpData::eIpStaticDataType::LOCAL_EQ_STRAIN:
            it.second.resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[eOutput::LOCAL_EQ_STRAIN] =
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::LOCAL_EQ_STRAIN);
            break;
        default:
            throw Exception(std::string("[") + __PRETTY_FUNCTION__ +
                                     "] this ip data type is not implemented.");
        }
    }
}


template <int TDim>
const Eigen::Vector3d NuTo::ContinuumBoundaryElement<TDim>::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    Eigen::VectorXd naturalSurfaceIpCoordinates = GetIntegrationType().GetLocalIntegrationPointCoordinates(rIpNum);

    Eigen::VectorXd naturalIpCoordinates =
            mInterpolationType->Get(Node::eDof::COORDINATES)
                    .CalculateNaturalSurfaceCoordinates(naturalSurfaceIpCoordinates, mSurfaceId);

    Eigen::VectorXd matrixN = mInterpolationType->Get(Node::eDof::COORDINATES).MatrixN(naturalIpCoordinates);
    Eigen::VectorXd nodeCoordinates = ExtractNodeValues(0, Node::eDof::COORDINATES);

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();
    globalIntegrationPointCoordinates.segment(0, GetLocalDimension()) = matrixN * nodeCoordinates;

    return globalIntegrationPointCoordinates;
}


#ifdef ENABLE_VISUALIZE
template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::Visualize(Visualize::UnstructuredGrid&,
                                                     const std::vector<eVisualizeWhat>&)
{
    std::cout << __PRETTY_FUNCTION__ << "Pleeeaaase, implement the visualization for me!!!" << std::endl;
}
#endif


namespace NuTo
{

template <int TDim>
int ContinuumBoundaryElement<TDim>::GetLocalDimension() const
{
    return TDim - 1;
}

template <int TDim>
int ContinuumBoundaryElement<TDim>::GetNumNodes() const
{
    return mInterpolationType->GetNumSurfaceNodes(mSurfaceId);
}

template <int TDim>
NodeBase* ContinuumBoundaryElement<TDim>::GetNode(int rLocalNodeNumber)
{
    int nodeId = mInterpolationType->GetSurfaceNodeIndex(mSurfaceId, rLocalNodeNumber);
    return const_cast<NodeBase*>(mBaseElement.GetNode(nodeId));
}

template <int TDim>
const NodeBase* ContinuumBoundaryElement<TDim>::GetNode(int rLocalNodeNumber) const
{
    int nodeId = mInterpolationType->GetSurfaceNodeIndex(mSurfaceId, rLocalNodeNumber);
    return mBaseElement.GetNode(nodeId);
}

template <int TDim>
int ContinuumBoundaryElement<TDim>::GetNumInfluenceNodes() const
{
    return mBaseElement.GetNumNodes();
}

template <int TDim>
const NodeBase* ContinuumBoundaryElement<TDim>::GetInfluenceNode(int rLocalNodeNumber) const
{
    return mBaseElement.GetNode(rLocalNodeNumber);
}

template <int TDim>
int ContinuumBoundaryElement<TDim>::GetNumNodes(Node::eDof rDofType) const
{
    return mInterpolationType->Get(rDofType).GetNumSurfaceNodes(mSurfaceId);
}

template <int TDim>
NodeBase* ContinuumBoundaryElement<TDim>::GetNode(int rLocalNodeNumber, Node::eDof rDofType)
{
    int nodeId = mInterpolationType->Get(rDofType).GetSurfaceNodeIndex(mSurfaceId, rLocalNodeNumber);
    return const_cast<NodeBase*>(mBaseElement.GetNode(nodeId));
}

template <int TDim>
const NodeBase* ContinuumBoundaryElement<TDim>::GetNode(int rLocalNodeNumber, Node::eDof rDofType) const
{
    int nodeId = mInterpolationType->Get(rDofType).GetSurfaceNodeIndex(mSurfaceId, rLocalNodeNumber);
    return mBaseElement.GetNode(nodeId);
}

template <int TDim>
std::shared_ptr<const Section> ContinuumBoundaryElement<TDim>::GetSection() const
{
    return mBaseElement.GetSection();
}

template <int TDim>
Eigen::VectorXd ContinuumBoundaryElement<TDim>::ExtractNodeValues(int rTimeDerivative, Node::eDof rDof) const
{
    return mBaseElement.ExtractNodeValues(rTimeDerivative, rDof);
}

template <int TDim>
double ContinuumBoundaryElement<TDim>::CalculateAlpha() const
{
    unsigned int theIP = 0; // This is a bit of a hack... I am sorry.

    double c = GetConstitutiveLaw(theIP).GetParameterDouble(Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS);
    return std::sqrt(c);
}


template <>
Eigen::Matrix<double, 0, 1> ContinuumBoundaryElement<1>::CalculateIPCoordinatesSurface(int) const
{
    return Eigen::Matrix<double, 0, 1>();
}

template <>
Eigen::Matrix<double, 1, 1> ContinuumBoundaryElement<2>::CalculateIPCoordinatesSurface(int rTheIP) const
{
    return GetIntegrationType().GetLocalIntegrationPointCoordinates(rTheIP);
}

template <>
Eigen::Matrix<double, 2, 1> ContinuumBoundaryElement<3>::CalculateIPCoordinatesSurface(int rTheIP) const
{
    return GetIntegrationType().GetLocalIntegrationPointCoordinates(rTheIP);
}


template <>
double NuTo::ContinuumBoundaryElement<1>::CalculateDetJxWeightIPxSection(double, int) const
{

    return mBaseElement.mSection->GetArea();
}

template <>
double NuTo::ContinuumBoundaryElement<2>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * GetIntegrationType().GetIntegrationPointWeight(rTheIP) *
           mBaseElement.mSection->GetThickness();
}

template <>
double NuTo::ContinuumBoundaryElement<3>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * GetIntegrationType().GetIntegrationPointWeight(rTheIP);
}

} // namespace NuTo

template class NuTo::ContinuumBoundaryElement<1>;
template class NuTo::ContinuumBoundaryElement<2>;
template class NuTo::ContinuumBoundaryElement<3>;
