#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

#include "nuto/mechanics/interpolationtypes/InterpolationType.h"

#include "nuto/mechanics/elements/EvaluateDataContinuum.h"

#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"

template<int TDim>
NuTo::ContinuumElement<TDim>::ContinuumElement(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIpDataType, rInterpolationType),
        mNodes(rNodes),
        mSection(nullptr)
{}

template<int TDim>
NuTo::Error::eError NuTo::ContinuumElement<TDim>::Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput)
{
    if (mSection == nullptr)
        throw MechanicsException(__PRETTY_FUNCTION__, "no section allocated for element.");

    EvaluateDataContinuum<TDim> data;
    ExtractAllNecessaryDofValues(data);

    auto constitutiveOutput = GetConstitutiveOutputMap(rElementOutput, data);

    auto constitutiveInput = GetConstitutiveInputMap(constitutiveOutput, data);
    constitutiveInput.insert(rInput.begin(), rInput.end());

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        CalculateNMatrixBMatrixDetJacobian(data, theIP);
        CalculateConstitutiveInputs(constitutiveInput, data);


        Error::eError error = EvaluateConstitutiveLaw<TDim>(constitutiveInput,constitutiveOutput,theIP);
        if (error != Error::SUCCESSFUL)
            return error;
        CalculateElementOutputs(rElementOutput, data, theIP);
    }
    return Error::SUCCESSFUL;
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::ExtractAllNecessaryDofValues(EvaluateDataContinuum<TDim>& data)
{
    // needs optimization,
    // not all dofs might be needed...

    const std::set<Node::eDof>& dofs = mInterpolationType->GetDofs();
    for (auto dof : dofs)
        if (mInterpolationType->IsConstitutiveInput(dof))
            data.mNodalValues[dof] = ExtractNodeValues(0, dof);

    data.mNodalValues[Node::COORDINATES] = ExtractNodeValues(0, Node::COORDINATES);

    if (mStructure->GetNumTimeDerivatives() >= 1)
        for (auto dof : dofs)
            if (mInterpolationType->IsConstitutiveInput(dof))
                data.mNodalValues_dt1[dof] = ExtractNodeValues(1, dof);
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElement<TDim>::ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);

    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofsPerNode = interpolationTypeDof.GetNumDofsPerNode();

    Eigen::VectorXd nodalValues(numDofsPerNode * numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase& node = *GetNode(iNode, rDofType);

        switch (rDofType)
        {
        case Node::COORDINATES:
            if (TDim == 1)
                nodalValues.segment<1>(iNode * 1) = node.GetCoordinates1D();
            if (TDim == 2)
                nodalValues.segment<2>(iNode * 2) = node.GetCoordinates2D();
            if (TDim == 3)
                nodalValues.segment<3>(iNode * 3) = node.GetCoordinates3D();
            break;
        case Node::DISPLACEMENTS:
            if (TDim == 1)
                nodalValues.segment<1>(iNode * 1) = node.GetDisplacements1D(rTimeDerivative);
            if (TDim == 2)
                nodalValues.segment<2>(iNode * 2) = node.GetDisplacements2D(rTimeDerivative);
            if (TDim == 3)
                nodalValues.segment<3>(iNode * 3) = node.GetDisplacements3D(rTimeDerivative);
            break;

        case Node::TEMPERATURE:
            nodalValues[iNode] = node.GetTemperature(rTimeDerivative);
            break;

        case Node::NONLOCALEQSTRAIN:
            nodalValues[iNode] = node.GetNonlocalEqStrain(rTimeDerivative);
            break;

        case Node::RELATIVEHUMIDITY:
            nodalValues[iNode] = node.GetRelativeHumidity(rTimeDerivative);
            break;

        case Node::WATERVOLUMEFRACTION:
            nodalValues[iNode] = node.GetWaterVolumeFraction(rTimeDerivative);
            break;

        case Node::DAMAGE:
            nodalValues[iNode] = node.GetDamage(rTimeDerivative);
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for " + Node::DofToString(rDofType));
        }
    }
    return nodalValues;
}

template<int TDim>
NuTo::ConstitutiveInputMap NuTo::ContinuumElement<TDim>::GetConstitutiveInputMap(
        const ConstitutiveOutputMap& rConstitutiveOutput,
        EvaluateDataContinuum<TDim>& rData) const
{
    ConstitutiveInputMap constitutiveInputMap = GetConstitutiveLaw(0)->GetConstitutiveInputs(rConstitutiveOutput, *GetInterpolationType());

    for (auto& itInput : constitutiveInputMap)
    {
        switch (itInput.first)
        {
        case Constitutive::Input::ENGINEERING_STRAIN:
            itInput.second = &(rData.mEngineeringStrain);
            break;

        case Constitutive::Input::NONLOCAL_EQ_STRAIN:
            itInput.second = &(rData.mNonlocalEqStrain);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY:
            itInput.second = &(rData.mRelativeHumidity);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY_DT1:
            itInput.second = &(rData.mRelativeHumidity_dt1);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT:
            itInput.second = &(rData.mRelativeHumidity_Gradient);
            break;

        case Constitutive::Input::TEMPERATURE:
            itInput.second = &(rData.mTemperature);
            break;

        case Constitutive::Input::TEMPERATURE_GRADIENT:
            itInput.second = &(rData.mTemperatureGradient);
            break;

        case Constitutive::Input::TEMPERATURE_CHANGE:
            itInput.second = &(rData.mTemperatureChange);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION:
            itInput.second = &(rData.mWaterVolumeFraction);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION_DT1:
            itInput.second = &(rData.mWaterVolumeFraction_dt1);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT:
            itInput.second = &(rData.mWaterVolumeFraction_Gradient);
            break;

        case Constitutive::Input::DAMAGE:
            itInput.second = &(rData.mDamage);
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive input " + Constitutive::InputToString(itInput.first) + " cannot be calculated by this element type.");
        }
    }
    return constitutiveInputMap;
}

template<int TDim>
NuTo::ConstitutiveOutputMap NuTo::ContinuumElement<TDim>::GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput, EvaluateDataContinuum<TDim> &rData) const
{
    ConstitutiveOutputMap constitutiveOutput;

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::INTERNAL_GRADIENT:
            FillConstitutiveOutputMapInternalGradient(constitutiveOutput, it.second->GetBlockFullVectorDouble(), rData);
            break;

        case Element::HESSIAN_0_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian0(constitutiveOutput, it.second->GetBlockFullMatrixDouble(), rData);
            break;

        case Element::HESSIAN_1_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian1(constitutiveOutput, it.second->GetBlockFullMatrixDouble(), rData);
            break;

        case Element::HESSIAN_2_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian2(constitutiveOutput, it.second->GetBlockFullMatrixDouble(), rData);
            break;

        case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
        {
            auto activeDofs = mInterpolationType->GetActiveDofs();
            if (activeDofs.size() > 1 && activeDofs.find(Node::DISPLACEMENTS) == activeDofs.end())
                throw MechanicsException(__PRETTY_FUNCTION__, "Lumped Hessian2 is only implemented for displacements.");
            int numDofs = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumDofs();
            it.second->GetBlockFullVectorDouble()[Node::DISPLACEMENTS].Resize(numDofs);
            rData.mTotalMass = 0;
            break;
        }

        case Element::UPDATE_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = nullptr;
            break;

        case Element::UPDATE_TMP_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA] = nullptr;
            break;

        case Element::IP_DATA:
            FillConstitutiveOutputMapIpData(constitutiveOutput, it.second->GetIpData(), rData);
            break;

        case Element::GLOBAL_ROW_DOF:
            CalculateGlobalRowDofs(it.second->GetBlockFullVectorInt());
            break;

        case Element::GLOBAL_COLUMN_DOF:
            CalculateGlobalColumnDofs(it.second->GetBlockFullVectorInt());
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
        }
    }
    return constitutiveOutput;
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapInternalGradient(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullVector<double>& rInternalGradient, EvaluateDataContinuum<TDim>& rData) const
{



    for (auto dofRow : mStructure->GetDofStatus().GetActiveDofTypes())
    {

        if (not (mInterpolationType->IsDof(dofRow)))
        {
            rInternalGradient[dofRow].Resize(0);
            continue;
        }

        rInternalGradient[dofRow].setZero(mInterpolationType->Get(dofRow).GetNumDofs());

        switch (dofRow)
        {
        case Node::DISPLACEMENTS:
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRESS] = &(rData.mEngineeringStress);
            break;
        case Node::NONLOCALEQSTRAIN:
            rConstitutiveOutput[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &(rData.mLocalEqStrain);
            rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &(rData.mNonlocalParameterXi);
            break;
        case Node::RELATIVEHUMIDITY:
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B] = &(rData.mInternalGradientRH_B);
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N] = &(rData.mInternalGradientRH_N);
            break;
        case Node::TEMPERATURE:
            rConstitutiveOutput[NuTo::Constitutive::Output::HEAT_FLUX] = &(rData.mHeatFlux);
            rConstitutiveOutput[NuTo::Constitutive::Output::HEAT_CHANGE] = &(rData.mHeatChange);
            break;
        case Node::WATERVOLUMEFRACTION:
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B] = &(rData.mInternalGradientWV_B);
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N] = &(rData.mInternalGradientWV_N);
            break;
        case Node::DAMAGE:
            rConstitutiveOutput[NuTo::Constitutive::Output::ELASTIC_ENERGY_DAMAGED_PART] = &(rData.mElasticEnergyDensity);
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");

        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapHessian0(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian0, EvaluateDataContinuum<TDim> &rData) const
{


    for (auto dofRow : mStructure->GetDofStatus().GetActiveDofTypes())
    {
        for (auto dofCol : mStructure->GetDofStatus().GetActiveDofTypes())
        {

            if (not (mInterpolationType->IsDof(dofRow) and mInterpolationType->IsDof(dofCol)))
            {
                rHessian0(dofRow, dofCol).Resize(0, 0);
                continue;
            }

            rHessian0(dofRow, dofCol).setZero(mInterpolationType->Get(dofRow).GetNumDofs(), mInterpolationType->Get(dofCol).GetNumDofs());

            if (not GetConstitutiveLaw(0)->CheckDofCombinationComputable(dofRow, dofCol, 0))
                continue;

            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::DISPLACEMENTS, Node::DISPLACEMENTS):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN] = &rData.mTangentStressStrain;
                break;
            case Node::CombineDofs(Node::DISPLACEMENTS, Node::NONLOCALEQSTRAIN):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN] = &rData.mTangentStressNonlocalEqStrain;
                break;
            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::DISPLACEMENTS):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN] = &rData.mTangentLocalEqStrainStrain;
                rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &rData.mNonlocalParameterXi;
                break;
            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::NONLOCALEQSTRAIN):
                rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &rData.mNonlocalParameterXi;
                break;
            case Node::CombineDofs(Node::TEMPERATURE, Node::TEMPERATURE):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT] = &rData.mTangentHeatFluxTemperatureGradient;
                break;

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY] = &rData.mEngineeringStress_dRH;
                break;

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::TEMPERATURE):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_TEMPERATURE] = &rData.mDStressDTemperature;
                break;

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION] = &rData.mEngineeringStress_dWV;
                break;

            case Node::CombineDofs(Node::RELATIVEHUMIDITY, Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0] = &rData.mInternalGradientRH_dRH_BB_H0;
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0] = &rData.mInternalGradientRH_dRH_NN_H0;
                break;

            case Node::CombineDofs(Node::RELATIVEHUMIDITY, Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0] = &rData.mInternalGradientRH_dWV_BN_H0;
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0] = &rData.mInternalGradientRH_dWV_NN_H0;
                break;

            case Node::CombineDofs(Node::WATERVOLUMEFRACTION, Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0] = &rData.mInternalGradientWV_dRH_NN_H0;
                break;

            case Node::CombineDofs(Node::WATERVOLUMEFRACTION, Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0] = &rData.mInternalGradientWV_dWV_BB_H0;
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0] = &rData.mInternalGradientWV_dWV_BN_H0;
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0] = &rData.mInternalGradientWV_dWV_NN_H0;
                break;

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::DAMAGE):
                rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRESS_DAMAGED_PART] = &rData.mEngineeringStressDamagedPart;
                break;

            case Node::CombineDofs(Node::DAMAGE, Node::DAMAGE):
                rConstitutiveOutput[NuTo::Constitutive::Output::ELASTIC_ENERGY_DAMAGED_PART] = &rData.mElasticEnergyDensity;
                break;

            case Node::CombineDofs(Node::DAMAGE, Node::DISPLACEMENTS):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN] = &(rData.mTangentElasticEnergyStrain);
                break;

            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_0_TIME_DERIVATIVE for (" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian1, EvaluateDataContinuum<TDim> &rData) const
{
    for (auto dofRow : mStructure->GetDofStatus().GetActiveDofTypes())
    {
        for (auto dofCol : mStructure->GetDofStatus().GetActiveDofTypes())
        {

            if (not (mInterpolationType->IsDof(dofRow) and mInterpolationType->IsDof(dofCol)))
            {
                rHessian1(dofRow, dofCol).Resize(0, 0);
                continue;
            }

            rHessian1(dofRow, dofCol).setZero(mInterpolationType->Get(dofRow).GetNumDofs(), mInterpolationType->Get(dofCol).GetNumDofs());

            if(!GetConstitutiveLaw(0)->CheckDofCombinationComputable(dofRow,dofCol,1))
                continue;

            switch (Node::CombineDofs(dofRow, dofCol))
            {

            case Node::CombineDofs(Node::RELATIVEHUMIDITY,Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1] = &rData.mInternalGradientRH_dRH_NN_H1;
                break;

            case Node::CombineDofs(Node::RELATIVEHUMIDITY,Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1] = &rData.mInternalGradientRH_dWV_NN_H1;
                break;

            case Node::CombineDofs(Node::WATERVOLUMEFRACTION,Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1] = &rData.mInternalGradientWV_dWV_NN_H1;
                break;

            case Node::CombineDofs(Node::TEMPERATURE, Node::TEMPERATURE):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_HEAT_D_TEMPERATURE] = &rData.mTangentHeatTemperature;
                break;
            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE):
            case Node::CombineDofs(Node::eDof::DAMAGE, Node::eDof::DAMAGE):
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_1_TIME_DERIVATIVE for (" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }

        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian2, EvaluateDataContinuum<TDim> &rData) const
{
    for (auto dofRow : mStructure->GetDofStatus().GetActiveDofTypes())
    {
        for (auto dofCol : mStructure->GetDofStatus().GetActiveDofTypes())
        {

            if (not (mInterpolationType->IsDof(dofRow) and mInterpolationType->IsDof(dofCol)))
            {
                rHessian2(dofRow, dofCol).Resize(0, 0);
                continue;
            }

            rHessian2(dofRow, dofCol).setZero(mInterpolationType->Get(dofRow).GetNumDofs(), mInterpolationType->Get(dofCol).GetNumDofs());

            if(!GetConstitutiveLaw(0)->CheckDofCombinationComputable(dofRow,dofCol,2))
                continue;

            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_2_TIME_DERIVATIVE for (" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapIpData(ConstitutiveOutputMap& rConstitutiveOutput, ElementOutputIpData& rIpData, EvaluateDataContinuum<TDim> &rData) const
{

    for (auto& it : rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    {                                       // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::DAMAGE:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::DAMAGE] = &(rData.mDamage);
            break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_VISUALIZE] = &(rData.mEngineeringPlasticStrainVisualize);
            break;
        case NuTo::IpData::ENGINEERING_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE] = &(rData.mEngineeringStrainVisualize);
            break;
        case NuTo::IpData::ENGINEERING_STRESS:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE] = &(rData.mEngineeringStressVisualize);
            break;
        case NuTo::IpData::EXTRAPOLATION_ERROR:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::EXTRAPOLATION_ERROR] = &(rData.mExtrapolationError);
            break;
        case NuTo::IpData::LOCAL_EQ_STRAIN:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &(rData.mLocalEqStrain);
            break;
        case NuTo::IpData::SHRINKAGE_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::SHRINKAGE_STRAIN_VISUALIZE] = &(rData.mShrinkageStrainVisualize);
            break;
        case NuTo::IpData::THERMAL_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::THERMAL_STRAIN] = &(rData.mThermalStrain);
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "this ip data type is not implemented.");
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateGlobalRowDofs(BlockFullVector<int> &rGlobalRowDofs) const
{
    const unsigned globalDimension = GetStructure()->GetDimension();

    for (auto dof : mStructure->GetDofStatus().GetActiveDofTypes())
    {

        if (not (mInterpolationType->IsDof(dof)))
        {
            rGlobalRowDofs[dof].Resize(0);
            continue;
        }

        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        const int numNodes = interpolationType.GetNumNodes();

        FullVector<int, Eigen::Dynamic>& dofWiseGlobalRowDofs = rGlobalRowDofs[dof];
        dofWiseGlobalRowDofs.setZero(interpolationType.GetNumDofs());

        switch (dof)
        {
        case Node::DISPLACEMENTS:
        {
            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];
                for (unsigned iDof = 0; iDof < globalDimension; ++iDof)
                    dofWiseGlobalRowDofs[globalDimension * iNodeDof + iDof] = nodePtr->GetDofDisplacement(iDof);
            }
            break;
        }
        case Node::TEMPERATURE:
        {
            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                dofWiseGlobalRowDofs[iNodeDof] = mNodes[interpolationType.GetNodeIndex(iNodeDof)]->GetDofTemperature();
            }
            break;
        }
        case Node::NONLOCALEQSTRAIN:
        {
            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                dofWiseGlobalRowDofs[iNodeDof] = mNodes[interpolationType.GetNodeIndex(iNodeDof)]->GetDofNonlocalEqStrain();
            }
            break;
        }
        case Node::RELATIVEHUMIDITY:
        {
            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                dofWiseGlobalRowDofs[iNodeDof] = mNodes[interpolationType.GetNodeIndex(iNodeDof)]->GetDofRelativeHumidity();
            }
            break;
        }
        case Node::WATERVOLUMEFRACTION:
        {
            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                dofWiseGlobalRowDofs[iNodeDof] = mNodes[interpolationType.GetNodeIndex(iNodeDof)]->GetDofWaterVolumeFraction();
            }
            break;
        }
        case Node::DAMAGE:
        {
            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                dofWiseGlobalRowDofs[iNodeDof] = mNodes[interpolationType.GetNodeIndex(iNodeDof)]->GetDofDamage();
            }
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for " + Node::DofToString(dof) + ".");
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateGlobalColumnDofs(BlockFullVector<int> &rGlobalDofMapping) const
{
    if (GetNumNonlocalElements() == 0)
        CalculateGlobalRowDofs(rGlobalDofMapping);
    else
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for nonlocal elements.");
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateConstitutiveInputs(const ConstitutiveInputMap& rConstitutiveInput, EvaluateDataContinuum<TDim> &rData)
{
    for (auto it : rConstitutiveInput)
    {
        switch (it.first)
        {
        case Constitutive::Input::ENGINEERING_STRAIN:
            rData.mEngineeringStrain.AsVector() = rData.mB.at(Node::DISPLACEMENTS) * rData.mNodalValues.at(Node::DISPLACEMENTS);
            break;

        case Constitutive::Input::NONLOCAL_EQ_STRAIN:
            rData.mNonlocalEqStrain.AsScalar() = (*rData.mN.at(Node::NONLOCALEQSTRAIN)) * rData.mNodalValues.at(Node::NONLOCALEQSTRAIN);
            break;

        case Constitutive::Input::DAMAGE:
            rData.mDamage.AsScalar() = (*rData.mN.at(Node::DAMAGE)) * rData.mNodalValues.at(Node::DAMAGE);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY:
            rData.mRelativeHumidity.AsScalar() = (*rData.mN.at(Node::RELATIVEHUMIDITY)) * rData.mNodalValues.at(Node::RELATIVEHUMIDITY);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY_DT1:
            rData.mRelativeHumidity_dt1.AsScalar() = (*rData.mN.at(Node::RELATIVEHUMIDITY)) * rData.mNodalValues_dt1.at(Node::RELATIVEHUMIDITY);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT:
            rData.mRelativeHumidity_Gradient.AsVector() = rData.mB.at(Node::RELATIVEHUMIDITY) * rData.mNodalValues.at(Node::RELATIVEHUMIDITY);
            break;

        case Constitutive::Input::TEMPERATURE:
            rData.mTemperature.AsScalar() = (*rData.mN.at(Node::TEMPERATURE)) * rData.mNodalValues.at(Node::TEMPERATURE);
            break;

        case Constitutive::Input::TEMPERATURE_GRADIENT:
            rData.mTemperatureGradient.AsVector() = rData.mB.at(Node::TEMPERATURE) * rData.mNodalValues.at(Node::TEMPERATURE);
            break;

        case Constitutive::Input::TEMPERATURE_CHANGE:
            if (mStructure->GetNumTimeDerivatives() >= 1)
                rData.mTemperatureChange.AsScalar() = (*rData.mN.at(Node::TEMPERATURE)) * rData.mNodalValues_dt1.at(Node::TEMPERATURE);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION:
            rData.mWaterVolumeFraction.AsScalar() = (*rData.mN.at(Node::WATERVOLUMEFRACTION)) * rData.mNodalValues.at(Node::WATERVOLUMEFRACTION);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION_DT1:
            rData.mWaterVolumeFraction_dt1.AsScalar() = (*rData.mN.at(Node::WATERVOLUMEFRACTION)) * rData.mNodalValues_dt1.at(Node::WATERVOLUMEFRACTION);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT:
            rData.mWaterVolumeFraction_Gradient.AsVector() = rData.mB.at(Node::WATERVOLUMEFRACTION) * rData.mNodalValues.at(Node::WATERVOLUMEFRACTION);
            break;

        case Constitutive::Input::TIME_STEP:
        case Constitutive::Input::CALCULATE_STATIC_DATA:
            break;


        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive input for " + Constitutive::InputToString(it.first) + " not implemented.");
        }
    }
}

template<int TDim>
Eigen::Matrix<double, TDim, TDim> NuTo::ContinuumElement<TDim>::CalculateJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::VectorXd& rNodeCoordinates) const
{
    int numCoordinateNodes = GetNumNodes(Node::COORDINATES);
    assert(rDerivativeShapeFunctions.rows() == numCoordinateNodes);
    assert(rDerivativeShapeFunctions.cols() == TDim);

    assert(rNodeCoordinates.rows() == TDim * GetNumNodes(Node::COORDINATES));

    Eigen::Matrix<double, TDim, Eigen::Dynamic> nodeBlockCoordinates(TDim, numCoordinateNodes);
    // convert the coordinates to a block structure
    // x0  x1  x1  x2 ...
    // y0  y1  y2  y3 ...
    // z0  z1  z2  z3 ...
    for (int i = 0; i < numCoordinateNodes; ++i)
        nodeBlockCoordinates.col(i) = rNodeCoordinates.block<TDim, 1>(TDim * i, 0);

    return nodeBlockCoordinates.lazyProduct(rDerivativeShapeFunctions);

}

template<int TDim>
Eigen::MatrixXd NuTo::ContinuumElement<TDim>::CalculateMatrixB(Node::eDof rDofType, const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::Matrix<double, TDim, TDim> rInvJacobian) const
{
    assert (rDerivativeShapeFunctions.rows() == GetNumNodes(rDofType));
    assert (rDerivativeShapeFunctions.cols() == TDim);


    Eigen::MatrixXd Bmat;
    switch (rDofType)
    {
    case Node::COORDINATES: // makes no sense, but is an example for gradient operator of vector valued dof type.
    {
        auto tmp = rDerivativeShapeFunctions.lazyProduct(rInvJacobian);
        /*   N0,x  N0,y  N0,z
         *   N1,x  N1,y  N1,z
         *   ...   ...   ...   */

        int numRows = rDerivativeShapeFunctions.rows();

        assert (tmp.cols() == TDim);
        assert (tmp.rows() == numRows);

        Bmat.resize(TDim, numRows*TDim);
        /*
         * transform to:
         *  N0,x   0    0    N1,x   0    0  ...
         *    0  N0,y   0      0  N1,y   0  ...
         *    0    0  N0,z     0    0  N1,z ...    */

        for (int i = 0; i < numRows; ++i)
            for (int iDim = 0; iDim < TDim; ++iDim)
                Bmat(iDim, TDim * i) = tmp(i, iDim);

        break;
    }
    case Node::DISPLACEMENTS:
    {
        Bmat = rDerivativeShapeFunctions.lazyProduct(rInvJacobian);
        BlowToBMatrixEngineeringStrain(Bmat);
        break;
    }
    default: // gradient for a scalar dof type
    {
        Bmat = rDerivativeShapeFunctions.lazyProduct(rInvJacobian).transpose();
        assert (Bmat.cols() == GetNumNodes(rDofType));
        assert (Bmat.rows() == TDim);
        break;
    }
    }

    return Bmat;
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput, EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    rData.mDetJxWeightIPxSection = CalculateDetJxWeightIPxSection(rData.mDetJacobian, rTheIP); // formerly known as "factor"

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::INTERNAL_GRADIENT:
            CalculateElementOutputInternalGradient(it.second->GetBlockFullVectorDouble(), rData, rTheIP);
            break;

        case Element::HESSIAN_0_TIME_DERIVATIVE:
            CalculateElementOutputHessian0(it.second->GetBlockFullMatrixDouble(), rData, rTheIP);
            break;

        case Element::HESSIAN_1_TIME_DERIVATIVE:
            CalculateElementOutputHessian1(it.second->GetBlockFullMatrixDouble(), rData, rTheIP);
            break;

        case Element::HESSIAN_2_TIME_DERIVATIVE:
            CalculateElementOutputHessian2(it.second->GetBlockFullMatrixDouble(), rData, rTheIP);
            break;

        case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
            for (auto dof : mInterpolationType->GetActiveDofs())
            {
                switch (dof)
                {
                case Node::DISPLACEMENTS:
                {
                    // calculate local mass matrix (the nonlocal terms are zero)
                    // don't forget to include determinant of the Jacobian and area
                    // detJ * area * density * HtH, :
                    FullVector<double, Eigen::Dynamic>& result = it.second->GetBlockFullVectorDouble()[Node::DISPLACEMENTS];
                    double rho = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY);
                    rData.mTotalMass += rData.mDetJxWeightIPxSection * rho;
                    const Eigen::VectorXd& shapeFunctions = mInterpolationType->Get(Node::DISPLACEMENTS).GetShapeFunctions(rTheIP);

                    //calculate for the translational dofs the diagonal entries
                    for (int i = 0; i < shapeFunctions.rows(); i++)
                       result(TDim * i) += shapeFunctions[i] * shapeFunctions[i] * rData.mDetJxWeightIPxSection * rho;

                    if (rTheIP + 1 == GetNumIntegrationPoints())
                    {
                        //calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
                        double sumDiagonal = result.Sum();

                        //scale so that the sum of the diagonals represents the full mass
                        double scaleFactor = rData.mTotalMass / sumDiagonal;
                        result *= scaleFactor;

                        // "spread" the entries over the whole result vector
                        for (int i = 0; i < shapeFunctions.rows(); i++)
                            for (int iDim = 1; iDim < TDim; ++iDim)
                                result(TDim * i + iDim) = result(TDim*i);
                    }
                }
                    break;
                default:
                    throw MechanicsException(__PRETTY_FUNCTION__, "Element output LUMPED_HESSIAN_2_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");
                }
            }
            break;

        case Element::UPDATE_STATIC_DATA:
        case Element::UPDATE_TMP_STATIC_DATA:
            break;
        case Element::IP_DATA:
            CalculateElementOutputIpData(it.second->GetIpData(), rData, rTheIP);
            break;
        case Element::GLOBAL_ROW_DOF:
        case Element::GLOBAL_COLUMN_DOF:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
        }
    }

}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputInternalGradient(BlockFullVector<double>& rInternalGradient, EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::DISPLACEMENTS:
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection *  rData.mB.at(dofRow).transpose() * rData.mEngineeringStress;
            break;

        case Node::NONLOCALEQSTRAIN:
        {
            const auto& N = *(rData.mN.at(dofRow));
            const auto& B = rData.mB.at(dofRow);
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection *
                                      ( N.transpose() * (rData.mNonlocalEqStrain[0] - rData.mLocalEqStrain[0]) / rData.mNonlocalParameterXi[0] +
                                        B.transpose() * (B * rData.mNodalValues.at(Node::NONLOCALEQSTRAIN)));
            break;
        }
        case Node::RELATIVEHUMIDITY:
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * (rData.mB.at(dofRow).transpose()  * rData.mInternalGradientRH_B +
                                                                         rData.mN.at(dofRow)->transpose() * rData.mInternalGradientRH_N);
            break;

        case Node::TEMPERATURE:
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * (rData.mB.at(dofRow).transpose() * rData.mHeatFlux +
                                                                         rData.mN.at(dofRow)->transpose() * rData.mHeatChange);
            break;

        case Node::WATERVOLUMEFRACTION:
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * (rData.mB.at(dofRow).transpose()  * rData.mInternalGradientWV_B +
                                                                         rData.mN.at(dofRow)->transpose() * rData.mInternalGradientWV_N);
            break;

        case Node::DAMAGE:
        {
            const auto  G       = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::FRACTURE_ENERGY);
            const auto  l       = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER);
            const auto& N       = *(rData.mN.at(Node::DAMAGE));
            const auto& B       = rData.mB.at(Node::DAMAGE);
            const auto& d       = rData.mNodalValues.at(Node::DAMAGE);
            const auto& d_dt    = rData.mNodalValues_dt1.at(Node::DAMAGE);

            const auto& kappa   = rData.mElasticEnergyDensity[0];
            const auto  visco   = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::ARTIFICIAL_VISCOSITY);

            //TODO: replace sigma * eps with static data kappa
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * (     ( (G/l + 2.*kappa) * N.transpose() * N
                                                                          +      G * l *B.transpose() * B
                                                                              )* d
                                                                          -      2. * N.transpose() * kappa
                                                                          +      N.transpose() * N * d_dt*visco
                                                                         );
            break;
        }


        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Element output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");
        }
    }

}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0, EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            if(!GetConstitutiveLaw(rTheIP)->CheckDofCombinationComputable(dofRow,dofCol,0))
                continue;
            auto& hessian0 = rHessian0(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
                hessian0 += rData.mDetJxWeightIPxSection *  rData.mB.at(dofRow).transpose() * rData.mTangentStressStrain * rData.mB.at(dofRow);
                break;

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN):
                hessian0 += rData.mDetJxWeightIPxSection *  rData.mB.at(dofRow).transpose() * rData.mTangentStressNonlocalEqStrain * *(rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::DISPLACEMENTS):
                hessian0 -= rData.mDetJxWeightIPxSection *  rData.mN.at(dofRow)->transpose() * rData.mTangentLocalEqStrainStrain.transpose() * (rData.mB.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::NONLOCALEQSTRAIN):
            {
                const auto& N = *(rData.mN.at(dofRow));
                const auto& B = rData.mB.at(dofRow);
                hessian0 += rData.mDetJxWeightIPxSection * (N.transpose() * (1./rData.mNonlocalParameterXi[0]) * N + B.transpose() * B);
                break;
            }

            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::TEMPERATURE):
                hessian0 += rData.mDetJxWeightIPxSection *  rData.mB.at(dofRow).transpose() * rData.mTangentHeatFluxTemperatureGradient * rData.mB.at(dofRow);
                break;

                //VHIRTHAMTODO get references to shape functions ---> no double find for the same value
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
                hessian0 += rData.mDetJxWeightIPxSection * (rData.mB.at(dofRow).transpose()  * rData.mInternalGradientRH_dRH_BB_H0[0] * rData.mB.at(dofCol) +
                                                            rData.mN.at(dofRow)->transpose() * rData.mInternalGradientRH_dRH_NN_H0 * (*rData.mN.at(dofCol)));
                break;

            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
                hessian0 += rData.mDetJxWeightIPxSection * (rData.mB.at(dofRow).transpose()  * rData.mInternalGradientRH_dWV_BN_H0 * (*rData.mN.at(dofCol)) +
                                                            rData.mN.at(dofRow)->transpose() * rData.mInternalGradientRH_dWV_NN_H0 * (*rData.mN.at(dofCol)));
                break;

            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
                hessian0 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow)->transpose()  * rData.mInternalGradientWV_dRH_NN_H0 * (*rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
                hessian0 += rData.mDetJxWeightIPxSection * (rData.mB.at(dofRow).transpose()  * rData.mInternalGradientWV_dWV_BB_H0[0] * rData.mB.at(dofCol) +
                                                            rData.mB.at(dofRow).transpose()  * rData.mInternalGradientWV_dWV_BN_H0 * (*rData.mN.at(dofCol)) +
                                                            rData.mN.at(dofRow)->transpose() * rData.mInternalGradientWV_dWV_NN_H0 * (*rData.mN.at(dofCol)));
                break;

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::RELATIVEHUMIDITY):
                hessian0 += rData.mDetJxWeightIPxSection * rData.mB.at(dofRow).transpose()  * rData.mEngineeringStress_dRH * (*rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE):
                hessian0 += rData.mDetJxWeightIPxSection * rData.mB.at(dofRow).transpose()  * rData.mDStressDTemperature * (*rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::WATERVOLUMEFRACTION):
                hessian0 += rData.mDetJxWeightIPxSection * rData.mB.at(dofRow).transpose()  * rData.mEngineeringStress_dWV * (*rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DAMAGE):
            {

                const auto& Bu = rData.mB.at(Node::eDof::DISPLACEMENTS);

                const auto& Nd = *(rData.mN.at(Node::eDof::DAMAGE));
                const double d = (Nd * rData.mNodalValues.at(Node::eDof::DAMAGE)).at(0,0);

                //TODO: substitute engineeringStress with tensile component
                hessian0 += rData.mDetJxWeightIPxSection * 2 * (d -1.) * Bu.transpose() * rData.mEngineeringStressDamagedPart * Nd;
                break;
            }

            case Node::CombineDofs(Node::eDof::DAMAGE, Node::eDof::DISPLACEMENTS):
            {
                const auto& Nd = *(rData.mN.at(Node::eDof::DAMAGE));
                const double d = (Nd * rData.mNodalValues.at(Node::eDof::DAMAGE)).at(0,0);

                const auto& Bu = rData.mB.at(Node::eDof::DISPLACEMENTS);

                //TODO: substitute engineeringStress with static data derivative d_kappa_d_strain
                hessian0 += rData.mDetJxWeightIPxSection * 2 * (d -1.) * Nd.transpose() * rData.mTangentElasticEnergyStrain.transpose() * Bu;
                break;
            }

            case Node::CombineDofs(Node::eDof::DAMAGE, Node::eDof::DAMAGE):
            {
                const auto  G = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::FRACTURE_ENERGY);
                const auto  l = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER);

                const auto& N = *(rData.mN.at(Node::eDof::DAMAGE));
                const auto& B = rData.mB.at(Node::eDof::DAMAGE);

                const auto& kappa   = rData.mElasticEnergyDensity[0];

                //TODO: replace sigma * eps with static data kappa
                hessian0 += rData.mDetJxWeightIPxSection * ( (G/l + 2*kappa) * N.transpose() * N +  G*l*B.transpose() * B );
                break;
            }

            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
//            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::DISPLACEMENTS):
//            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::DISPLACEMENTS):
//            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::DISPLACEMENTS):
//                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for (" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputHessian1(BlockFullMatrix<double>& rHessian1, EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            if(!GetConstitutiveLaw(rTheIP)->CheckDofCombinationComputable(dofRow,dofCol,1))
                continue;
            auto& hessian1 = rHessian1(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
                break;

            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::TEMPERATURE):
                hessian1 += rData.mDetJxWeightIPxSection
                          * rData.mN.at(dofRow)->transpose()
                          * rData.mTangentHeatTemperature
                          * (*rData.mN.at(dofCol));
                break;
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
                hessian1 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow)->transpose() * rData.mInternalGradientRH_dRH_NN_H1 * (*rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
                hessian1 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow)->transpose() * rData.mInternalGradientRH_dWV_NN_H1 * (*rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
                break;

            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
                hessian1 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow)->transpose() * rData.mInternalGradientWV_dWV_NN_H1 * (*rData.mN.at(dofCol));
                break;

            case Node::CombineDofs(Node::eDof::DAMAGE, Node::eDof::DAMAGE):
            {
                const auto& N = *(rData.mN.at(Node::eDof::DAMAGE));
                const auto  visco   = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::ARTIFICIAL_VISCOSITY);
                hessian1 += visco * rData.mDetJxWeightIPxSection * N.transpose() * N;
            }
                break;

            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
//            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::RELATIVEHUMIDITY):
//            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::WATERVOLUMEFRACTION):
//            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::DISPLACEMENTS):
//            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::DISPLACEMENTS):
//            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::DISPLACEMENTS):
//            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE):
//                break;
            default:
                throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Element output HESSIAN_1_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputHessian2(BlockFullMatrix<double>& rHessian2, EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            if(!GetConstitutiveLaw(rTheIP)->CheckDofCombinationComputable(dofRow,dofCol,2))
                continue;
            auto& hessian2 = rHessian2(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            {
                const auto& N = *(rData.mN.at(dofRow));
                double rho = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY);
                hessian2 += rho * N.transpose() * N * rData.mDetJxWeightIPxSection;

                break;
            }
            default:
                throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Element output HESSIAN_2_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputIpData(ElementOutputIpData& rIpData, EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    for (auto& it : rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    {                                       // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::DAMAGE:
            it.second.col(rTheIP) = std::move(rData.mDamage);
            break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mEngineeringPlasticStrainVisualize);
            break;
        case NuTo::IpData::ENGINEERING_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mEngineeringStrainVisualize);
            break;
        case NuTo::IpData::ENGINEERING_STRESS:
            it.second.col(rTheIP) = std::move(rData.mEngineeringStressVisualize);
            break;
        case NuTo::IpData::EXTRAPOLATION_ERROR:
            it.second.col(rTheIP) = std::move(rData.mExtrapolationError);
            break;
        case NuTo::IpData::LOCAL_EQ_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mLocalEqStrain);
            break;
        case NuTo::IpData::SHRINKAGE_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mShrinkageStrainVisualize);
            break;
        case NuTo::IpData::THERMAL_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mThermalStrain);
            break;
        default:
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Ip data not implemented.");
        }
    }
}

template<int TDim>
NuTo::NodeBase* NuTo::ContinuumElement<TDim>::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

template<int TDim>
const NuTo::NodeBase* NuTo::ContinuumElement<TDim>::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

template<int TDim>
NuTo::NodeBase* NuTo::ContinuumElement<TDim>::GetNode(int rLocalNodeNumber, Node::eDof rDofType)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

template<int TDim>
const NuTo::NodeBase* NuTo::ContinuumElement<TDim>::GetNode(int rLocalNodeNumber, Node::eDof rDofType) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    assert(mNodes[nodeIndex] != nullptr);
    return mNodes[nodeIndex];
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    assert(rNode != nullptr);
    mNodes[rLocalNodeNumber] = rNode;
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::ResizeNodes(int rNewNumNodes)
{
    if (rNewNumNodes == (int) mNodes.size())
        return;

    if (rNewNumNodes > (int) mNodes.size())
    {
        // just resize (enlarge)
        mNodes.resize(rNewNumNodes);
    }
    else
    {
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Resize that reduces the number of nodes is not implemented yet.");
    }
}

template<int TDim>
const Eigen::VectorXd NuTo::ContinuumElement<TDim>::GetIntegrationPointVolume() const
{
    Eigen::MatrixXd nodeCoordinates = this->ExtractNodeValues(0, Node::COORDINATES);

    Eigen::VectorXd volume(GetNumIntegrationPoints());
    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        Eigen::MatrixXd derivativeShapeFunctionsNatural = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);
        double detJacobian = CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates).determinant();
        volume[theIP] = detJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
    return volume;
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CheckElement()
{
    int numIntegrationPoints = GetNumIntegrationPoints();

    if (numIntegrationPoints < 1)
    {
        MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid integration type.");
    }

    int theIP = 0;
    const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);
    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    double detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates).determinant();
    if (detJacobian < 0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    }

    double size = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
        detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates).determinant();
        if (detJacobian <= 0)
        {
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Determinant of the Jacobian <= zero, no inversion possible.");
        }
        size += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    assert(std::abs(GetIntegrationPointVolume().sum() / size - 1) < 1.e-6); // just to be sure ...

    // check element length
    if (size < 1e-14)
    {
        MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] element with zero size (check nodes).");
    }

}

template<int TDim>
void NuTo::ContinuumElement<TDim>::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    std::replace(mNodes.begin(), mNodes.end(), rOldPtr, rNewPtr);
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateNMatrixBMatrixDetJacobian(EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    // calculate Jacobian
    const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(rTheIP);

    Eigen::Matrix<double, TDim, TDim> jacobian = CalculateJacobian(derivativeShapeFunctionsGeometryNatural, rData.mNodalValues[Node::COORDINATES]);
    rData.mDetJacobian = jacobian.determinant();
    if (rData.mDetJacobian == 0)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Determinant of the Jacobian is zero, no inversion possible.");

    Eigen::Matrix<double, TDim, TDim> invJacobian = jacobian.inverse();

    // calculate shape functions and their derivatives
    for (auto dof : mInterpolationType->GetDofs())
    {
        if (dof == Node::COORDINATES)
            continue;
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        rData.mN[dof] = &interpolationType.GetMatrixN(rTheIP);

        rData.mB[dof] = CalculateMatrixB(dof, interpolationType.GetDerivativeShapeFunctionsNatural(rTheIP), invJacobian);
    }
}

namespace NuTo // template specialization in *.cpp somehow requires the definition to be in the namespace...
{
template<>
void NuTo::ContinuumElement<1>::BlowToBMatrixEngineeringStrain(Eigen::MatrixXd& rDerivativeShapeFunctions) const
{
    assert(rDerivativeShapeFunctions.cols() == 1);
    rDerivativeShapeFunctions.transposeInPlace();
}

template<>
void NuTo::ContinuumElement<2>::BlowToBMatrixEngineeringStrain(Eigen::MatrixXd& rDerivativeShapeFunctions) const
{
    int numNodes = GetNumNodes(Node::DISPLACEMENTS);
    assert(rDerivativeShapeFunctions.cols() == 2);
    assert(rDerivativeShapeFunctions.rows() == numNodes);

    Eigen::MatrixXd Bmat = Eigen::MatrixXd::Zero(3, numNodes * 2);

    for (int iNode = 0, iColumn = 0; iNode < numNodes; ++iNode, iColumn += 2)
    {
        double dNdX = rDerivativeShapeFunctions(iNode, 0);
        double dNdY = rDerivativeShapeFunctions(iNode, 1);

        Bmat(0, iColumn) = dNdX;
        Bmat(1, iColumn + 1) = dNdY;
        Bmat(2, iColumn) = dNdY;
        Bmat(2, iColumn + 1) = dNdX;
    }
    std::swap(rDerivativeShapeFunctions, Bmat);
}

template<>
void NuTo::ContinuumElement<3>::BlowToBMatrixEngineeringStrain(Eigen::MatrixXd& rDerivativeShapeFunctions) const
{
    int numNodes = GetNumNodes(Node::DISPLACEMENTS);
    assert(rDerivativeShapeFunctions.cols() == 3);
    assert(rDerivativeShapeFunctions.rows() == numNodes);
    Eigen::MatrixXd Bmat = Eigen::MatrixXd::Zero(6, numNodes * 3);

    for (int iNode = 0, iColumn = 0; iNode < numNodes; ++iNode, iColumn += 3)
    {
        double dNdX = rDerivativeShapeFunctions(iNode, 0);
        double dNdY = rDerivativeShapeFunctions(iNode, 1);
        double dNdZ = rDerivativeShapeFunctions(iNode, 2);

        /* according to Jirásek
         *
         *     +0  +1  +2
         *    -------------
         * 0 |  dx  0   0  |     - e_x
         * 1 |  0   dy  0  |     - e_y
         * 2 |  0   0   dz |     - e_z
         * 3 |  0   dz  dy |     - g_yz
         * 4 |  dz  0   dx |     - g_xz
         * 5 |  dy  dx  0  |     - g_xy
         *    -------------
         */


        Bmat(0, iColumn)     = dNdX;
        Bmat(1, iColumn + 1) = dNdY;
        Bmat(2, iColumn + 2) = dNdZ;

        Bmat(3, iColumn + 1) = dNdZ;
        Bmat(3, iColumn + 2) = dNdY;

        Bmat(4, iColumn)     = dNdZ;
        Bmat(4, iColumn + 2) = dNdX;

        Bmat(5, iColumn)     = dNdY;
        Bmat(5, iColumn + 1) = dNdX;
    }
    std::swap(rDerivativeShapeFunctions, Bmat);
}

template<>
double NuTo::ContinuumElement<1>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    Eigen::MatrixXd matrixN = mInterpolationType->Get(Node::COORDINATES).GetMatrixN(rTheIP);
    Eigen::VectorXd globalIPCoordinate = matrixN * ExtractNodeValues(0, Node::COORDINATES);

    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetArea() * mSection->AsSectionTruss()->GetAreaFactor(globalIPCoordinate.at(0, 0));
}

template<>
double NuTo::ContinuumElement<2>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetThickness();
}

template<>
double NuTo::ContinuumElement<3>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP);
}

template<>
NuTo::ConstitutiveStaticDataBase* ContinuumElement<1>::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticData1D(this);
}
template<>
NuTo::ConstitutiveStaticDataBase* ContinuumElement<2>::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticData2D(this);
}
template<>
NuTo::ConstitutiveStaticDataBase* ContinuumElement<3>::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticData3D(this);
}

template<>
const ContinuumElement<1>& ContinuumElement<1>::AsContinuumElement1D() const
{
    return *this;
}

template<>
const ContinuumElement<2>& ContinuumElement<2>::AsContinuumElement2D() const
{
    return *this;
}

template<>
const ContinuumElement<3>& ContinuumElement<3>::AsContinuumElement3D() const
{
    return *this;
}

template<>
ContinuumElement<1>& ContinuumElement<1>::AsContinuumElement1D()
{
    return *this;
}

template<>
ContinuumElement<2>& ContinuumElement<2>::AsContinuumElement2D()
{
    return *this;
}

template<>
ContinuumElement<3>& ContinuumElement<3>::AsContinuumElement3D()
{
    return *this;
}

}  // namespace NuTo

template class NuTo::ContinuumElement<1>;
template class NuTo::ContinuumElement<2>;
template class NuTo::ContinuumElement<3>;


#ifdef ENABLE_SERIALIZATION
template void NuTo::ContinuumElement<1>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<2>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<3>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<1>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<2>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<3>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<1>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<2>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<3>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElement<TDim>::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ContinuumElement " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElement_ElementBase",boost::serialization::base_object<ElementBase >(*this));
    ar & boost::serialization::make_nvp("mSection", const_cast<SectionBase*&>(mSection));

    const std::uintptr_t* mNodesAddress = reinterpret_cast<const std::uintptr_t*>(mNodes.data());
    int size = mNodes.size();
    ar & boost::serialization::make_nvp("mNodes_size", size);
    ar & boost::serialization::make_nvp("mNodes", boost::serialization::make_array(mNodesAddress, size));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ContinuumElement" << std::endl;
#endif
}

template void NuTo::ContinuumElement<1>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<2>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<3>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<1>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<2>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<3>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<1>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<2>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElement<3>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElement<TDim>::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start deserialize ContinuumElement " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElement_ElementBase",boost::serialization::base_object<ElementBase >(*this));
    ar & boost::serialization::make_nvp("mSection", const_cast<SectionBase*&>(mSection));

    int size = 0;
    ar & boost::serialization::make_nvp("mNodes_size", size);
    std::uintptr_t* mNodesAddress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mNodes", boost::serialization::make_array(mNodesAddress, size));
    mNodes.assign(reinterpret_cast<NodeBase**>(&mNodesAddress[0]), reinterpret_cast<NodeBase**>(&mNodesAddress[size]));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish deserialize ContinuumElement" << std::endl;
#endif
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
{
    for(std::vector<NodeBase*>::iterator it = mNodes.begin(); it != mNodes.end(); it++)
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(*it));
        if (itCast!=mNodeMapCast.end())
        {
            *it = reinterpret_cast<NodeBase*>(itCast->second);
        }
        else
            throw MechanicsException(__PRETTY_FUNCTION__, "The NodeBase-Pointer could not be updated.");
    }
}

BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElement<1>)), "ContinuumElement_1")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElement<2>)), "ContinuumElement_2")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElement<3>)), "ContinuumElement_3")
#endif // ENABLE_SERIALIZATION
