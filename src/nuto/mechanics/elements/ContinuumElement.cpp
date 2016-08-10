#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

#include "nuto/mechanics/interpolationtypes/InterpolationType.h"

#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

template<int TDim>
NuTo::ContinuumElement<TDim>::ContinuumElement(const NuTo::StructureBase* rStructure,
        const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType,
        IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIpDataType, rInterpolationType),
        mNodes(rNodes),
        mSection(nullptr)
{}

template<int TDim>
NuTo::Error::eError NuTo::ContinuumElement<TDim>::Evaluate(const ConstitutiveInputMap& rInput,
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput)
{
    if (mSection == nullptr)
        throw MechanicsException(__PRETTY_FUNCTION__, "no section allocated for element.");

    ExtractAllNecessaryDofValues();

    auto constitutiveOutput = GetConstitutiveOutputMap(rElementOutput);
    auto constitutiveInput = GetConstitutiveInputMap(constitutiveOutput);

    constitutiveInput.Merge(rInput);

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        CalculateNMatrixBMatrixDetJacobian(theIP);
        CalculateConstitutiveInputs(constitutiveInput);

        mDetJxWeightIPxSection = CalculateDetJxWeightIPxSection(theIP); // formerly known as "factor"
        Error::eError error = EvaluateConstitutiveLaw<TDim>(constitutiveInput,constitutiveOutput,theIP);
        if (error != Error::SUCCESSFUL)
            return error;
        CalculateElementOutputs(rElementOutput, theIP, constitutiveInput, constitutiveOutput);
    }
    return Error::SUCCESSFUL;
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::ExtractAllNecessaryDofValues()
{
    //! @todo needs optimization, not all dofs might be needed...

    const std::set<Node::eDof>& dofs = mInterpolationType->GetDofs();
    for (auto dof : dofs)
        if (mInterpolationType->IsConstitutiveInput(dof))
            mNodalValues[dof] = ExtractNodeValues(0, dof);

    mNodalValues[Node::COORDINATES] = ExtractNodeValues(0, Node::COORDINATES);

    if (mStructure->GetNumTimeDerivatives() >= 1)
        for (auto dof : dofs)
            if (mInterpolationType->IsConstitutiveInput(dof))
                mNodalValues_dt1[dof] = ExtractNodeValues(1, dof);
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

        nodalValues.block(iNode * numDofsPerNode, 0, numDofsPerNode, 1) = node.Get(rDofType, rTimeDerivative);

    }

    return nodalValues;
}

template<int TDim>
NuTo::ConstitutiveInputMap NuTo::ContinuumElement<TDim>::GetConstitutiveInputMap(
        const ConstitutiveOutputMap& rConstitutiveOutput) const
{
    // create maps with only the keys
    ConstitutiveInputMap constitutiveInputMap = GetConstitutiveLaw(0)->GetConstitutiveInputs(rConstitutiveOutput, *GetInterpolationType());

    // attach corresponding scalar/vector/matrix object to each key
    for (auto& itInput : constitutiveInputMap)
    {
        itInput.second = ConstitutiveIOBase::makeConstitutiveIO<TDim>(itInput.first);
    }
    return constitutiveInputMap;
}

template<int TDim>
NuTo::ConstitutiveOutputMap NuTo::ContinuumElement<TDim>::GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) const
{
    ConstitutiveOutputMap constitutiveOutput;

    // find the outputs we need
    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::INTERNAL_GRADIENT:
            FillConstitutiveOutputMapInternalGradient(constitutiveOutput, it.second->GetBlockFullVectorDouble());
            break;

        case Element::HESSIAN_0_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian0(constitutiveOutput, it.second->GetBlockFullMatrixDouble());
            break;

        case Element::HESSIAN_1_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian1(constitutiveOutput, it.second->GetBlockFullMatrixDouble());
            break;

        case Element::HESSIAN_2_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian2(constitutiveOutput, it.second->GetBlockFullMatrixDouble());
            break;

        case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
        {
            auto activeDofs = mInterpolationType->GetActiveDofs();
            if (activeDofs.size() > 1 && activeDofs.find(Node::DISPLACEMENTS) == activeDofs.end())
                throw MechanicsException(__PRETTY_FUNCTION__, "Lumped Hessian2 is only implemented for displacements.");
            int numDofs = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumDofs();
            it.second->GetBlockFullVectorDouble()[Node::DISPLACEMENTS].Resize(numDofs);
            mTotalMass = 0;
            break;
        }

        case Element::UPDATE_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = nullptr;
            break;

        case Element::UPDATE_TMP_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA] = nullptr;
            break;

        case Element::IP_DATA:
            FillConstitutiveOutputMapIpData(constitutiveOutput, it.second->GetIpData());
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

    // allocate the objects for the output data
    for (auto& outputs : constitutiveOutput)
    {
        outputs.second = ConstitutiveIOBase::makeConstitutiveIO<TDim>(outputs.first);
    }
    return constitutiveOutput;
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapInternalGradient(
        ConstitutiveOutputMap& rConstitutiveOutput,
        BlockFullVector<double>& rInternalGradient) const
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
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRESS];
            break;
        case Node::NONLOCALEQSTRAIN:
            rConstitutiveOutput[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN];
            rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI];
            break;
        case Node::RELATIVEHUMIDITY:
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B];
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N];
            break;
        case Node::TEMPERATURE:
            rConstitutiveOutput[NuTo::Constitutive::Output::HEAT_FLUX];
            rConstitutiveOutput[NuTo::Constitutive::Output::HEAT_CHANGE];
            break;
        case Node::WATERVOLUMEFRACTION:
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B];
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N];
            break;
        case Node::CRACKPHASEFIELD:
            rConstitutiveOutput[NuTo::Constitutive::Output::ELASTIC_ENERGY_DAMAGED_PART];
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");

        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapHessian0(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian0) const
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
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN];
                break;
            case Node::CombineDofs(Node::DISPLACEMENTS, Node::NONLOCALEQSTRAIN):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN];
                break;
            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::DISPLACEMENTS):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN];
                rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI];
                break;
            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::NONLOCALEQSTRAIN):
                rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI];
                break;
            case Node::CombineDofs(Node::TEMPERATURE, Node::TEMPERATURE):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT];
                break;

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY];
                break;

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::TEMPERATURE):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_TEMPERATURE];
                break;

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION];
                break;

            case Node::CombineDofs(Node::RELATIVEHUMIDITY, Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0];
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0];
                break;

            case Node::CombineDofs(Node::RELATIVEHUMIDITY, Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0];
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0];
                break;

            case Node::CombineDofs(Node::WATERVOLUMEFRACTION, Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0];
                break;

            case Node::CombineDofs(Node::WATERVOLUMEFRACTION, Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0];
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0];
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0];
                break;

            default:
                throw MechanicsException(__PRETTY_FUNCTION__,
                        "Constitutive output HESSIAN_0_TIME_DERIVATIVE for ("
                        + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian1) const
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
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1];
                break;

            case Node::CombineDofs(Node::RELATIVEHUMIDITY,Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1];
                break;

            case Node::CombineDofs(Node::WATERVOLUMEFRACTION,Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1];
                break;

            case Node::CombineDofs(Node::TEMPERATURE, Node::TEMPERATURE):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_HEAT_D_TEMPERATURE];
                break;
            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE):
            case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::CRACKPHASEFIELD):
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_1_TIME_DERIVATIVE for (" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }

        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian2) const
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
void NuTo::ContinuumElement<TDim>::FillConstitutiveOutputMapIpData(ConstitutiveOutputMap& rConstitutiveOutput, ElementOutputIpData& rIpData) const
{

    for (auto& it : rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    {                                       // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::DAMAGE:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::DAMAGE];
            break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_VISUALIZE];
            break;
        case NuTo::IpData::ENGINEERING_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE];
            break;
        case NuTo::IpData::ENGINEERING_STRESS:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE];
            break;
        case NuTo::IpData::EXTRAPOLATION_ERROR:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::EXTRAPOLATION_ERROR];
            break;
        case NuTo::IpData::LOCAL_EQ_STRAIN:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN];
            break;
        case NuTo::IpData::SHRINKAGE_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::SHRINKAGE_STRAIN_VISUALIZE];
            break;
        case NuTo::IpData::THERMAL_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::THERMAL_STRAIN];
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "this ip data type is not implemented.");
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateGlobalRowDofs(BlockFullVector<int> &rGlobalRowDofs) const
{
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


        unsigned int numDofsPerType = mNodes[interpolationType.GetNodeIndex(0)]->GetNum(dof);

        for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
        {
            const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];

            for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
            {
                dofWiseGlobalRowDofs[numDofsPerType * iNodeDof + iDof] = nodePtr->GetDof(dof, iDof);
            }
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
void NuTo::ContinuumElement<TDim>::CalculateConstitutiveInputs(ConstitutiveInputMap& rConstitutiveInput)
{
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    for (auto& it : rConstitutiveInput)
    {
        switch (it.first)
        {
        case Constitutive::Input::ENGINEERING_STRAIN:        
        {
            auto& strain = *static_cast<ConstitutiveVector<VoigtDim>*>(it.second.get());
            strain.AsVector() = mB.at(Node::DISPLACEMENTS) * mNodalValues.at(Node::DISPLACEMENTS);
            break;
        }
        case Constitutive::Input::NONLOCAL_EQ_STRAIN:
        {
            auto& nonLocalEqStrain = *static_cast<ConstitutiveScalar*>(it.second.get());
            nonLocalEqStrain.AsScalar() = (*mN.at(Node::NONLOCALEQSTRAIN)) * mNodalValues.at(Node::NONLOCALEQSTRAIN);
            break;
        }
        case Constitutive::Input::CRACK_PHASE_FIELD:
        {
            auto& damage = *static_cast<ConstitutiveScalar*>(it.second.get());
            damage.AsScalar() = (*mN.at(Node::CRACKPHASEFIELD)) * mNodalValues.at(Node::CRACKPHASEFIELD);
            break;
        }
        case Constitutive::Input::RELATIVE_HUMIDITY:
        {
            auto& relativeHumidity = *static_cast<ConstitutiveScalar*>(it.second.get());
            relativeHumidity.AsScalar() = (*mN.at(Node::RELATIVEHUMIDITY)) * mNodalValues.at(Node::RELATIVEHUMIDITY);
            break;
        }
        case Constitutive::Input::RELATIVE_HUMIDITY_DT1:
        {
            auto& relativeHumidity_dt1 = *static_cast<ConstitutiveScalar*>(it.second.get());
            relativeHumidity_dt1.AsScalar() = (*mN.at(Node::RELATIVEHUMIDITY)) * mNodalValues_dt1.at(Node::RELATIVEHUMIDITY);
            break;
        }
        case Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT:
        {
            auto& relHumidityGrad = *static_cast<ConstitutiveVector<TDim>*>(it.second.get());
            relHumidityGrad.AsVector() = mB.at(Node::RELATIVEHUMIDITY) * mNodalValues.at(Node::RELATIVEHUMIDITY);
            break;
        }
        case Constitutive::Input::TEMPERATURE:
        {
            auto& temperature = *static_cast<ConstitutiveScalar*>(it.second.get());
            temperature.AsScalar() = (*mN.at(Node::TEMPERATURE)) * mNodalValues.at(Node::TEMPERATURE);
            break;
        }
        case Constitutive::Input::TEMPERATURE_GRADIENT:
        {
            auto& tempGradient = *static_cast<ConstitutiveVector<TDim>*>(it.second.get());   
            tempGradient.AsVector() = mB.at(Node::TEMPERATURE) * mNodalValues.at(Node::TEMPERATURE);
            break;
        }
        case Constitutive::Input::TEMPERATURE_CHANGE:
        {
            if (mStructure->GetNumTimeDerivatives() >= 1)
            {
                auto& temperatureChange = *static_cast<ConstitutiveScalar*>(it.second.get());
                temperatureChange.AsScalar() = (*mN.at(Node::TEMPERATURE)) * mNodalValues_dt1.at(Node::TEMPERATURE);
            }
            break;
        }
        case Constitutive::Input::WATER_VOLUME_FRACTION:
        {
            auto& waterVolumeFraction = *static_cast<ConstitutiveScalar*>(it.second.get());
            waterVolumeFraction.AsScalar() = (*mN.at(Node::WATERVOLUMEFRACTION)) * mNodalValues.at(Node::WATERVOLUMEFRACTION);
            break;
        }
        case Constitutive::Input::WATER_VOLUME_FRACTION_DT1:
        {
            auto& waterVolumeFraction_dt1 = *static_cast<ConstitutiveScalar*>(it.second.get());
            waterVolumeFraction_dt1.AsScalar() = (*mN.at(Node::WATERVOLUMEFRACTION)) * mNodalValues_dt1.at(Node::WATERVOLUMEFRACTION);
            break;
        }
        case Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT:
        {
            auto& waterVolumeFractionGrad = *static_cast<ConstitutiveVector<TDim>*>(it.second.get());
            waterVolumeFractionGrad.AsVector() = mB.at(Node::WATERVOLUMEFRACTION) * mNodalValues.at(Node::WATERVOLUMEFRACTION);
            break;
        }
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
void NuTo::ContinuumElement<TDim>::CalculateElementOutputs(
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput, int rTheIP,
        const ConstitutiveInputMap& constitutiveInputMap,
        const ConstitutiveOutputMap& constitutiveOutputMap) const
{

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::INTERNAL_GRADIENT:
            CalculateElementOutputInternalGradient(it.second->GetBlockFullVectorDouble(), rTheIP, constitutiveInputMap, constitutiveOutputMap);
            break;

        case Element::HESSIAN_0_TIME_DERIVATIVE:
            CalculateElementOutputHessian0(it.second->GetBlockFullMatrixDouble(), rTheIP, constitutiveOutputMap);
            break;

        case Element::HESSIAN_1_TIME_DERIVATIVE:
            CalculateElementOutputHessian1(it.second->GetBlockFullMatrixDouble(), rTheIP, constitutiveOutputMap);
            break;

        case Element::HESSIAN_2_TIME_DERIVATIVE:
            CalculateElementOutputHessian2(it.second->GetBlockFullMatrixDouble(), rTheIP);
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
                    mTotalMass += mDetJxWeightIPxSection * rho;
                    const Eigen::VectorXd& shapeFunctions = mInterpolationType->Get(Node::DISPLACEMENTS).GetShapeFunctions(rTheIP);

                    //calculate for the translational dofs the diagonal entries
                    for (int i = 0; i < shapeFunctions.rows(); i++)
                       result(TDim * i) += shapeFunctions[i] * shapeFunctions[i] * mDetJxWeightIPxSection * rho;

                    if (rTheIP + 1 == GetNumIntegrationPoints())
                    {
                        //calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
                        double sumDiagonal = result.Sum();

                        //scale so that the sum of the diagonals represents the full mass
                        double scaleFactor = mTotalMass / sumDiagonal;
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
            CalculateElementOutputIpData(it.second->GetIpData(), rTheIP, constitutiveOutputMap);
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
void NuTo::ContinuumElement<TDim>::CalculateElementOutputInternalGradient(
        BlockFullVector<double>& rInternalGradient, int rTheIP,
        const ConstitutiveInputMap& constitutiveInputMap,
        const ConstitutiveOutputMap& constitutiveOutputMap) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::DISPLACEMENTS:
        {
            const auto& engineeringStress= *static_cast<EngineeringStress<TDim>*>(constitutiveOutputMap.at(Constitutive::Output::ENGINEERING_STRESS).get());
            rInternalGradient[dofRow] += mDetJxWeightIPxSection *  mB.at(dofRow).transpose() * engineeringStress;
            break;
        }
        case Node::NONLOCALEQSTRAIN:
        {
            const auto& N = *(mN.at(dofRow));
            const auto& B = mB.at(dofRow);
            const auto& nonlocalEqStrain = *static_cast<ConstitutiveScalar*>(constitutiveInputMap.at(Constitutive::Input::NONLOCAL_EQ_STRAIN).get());
            const auto& localEqStrain = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::LOCAL_EQ_STRAIN).get());
            const auto& nonlocalXi = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::NONLOCAL_PARAMETER_XI).get());
            rInternalGradient[dofRow] += mDetJxWeightIPxSection *
                                      ( N.transpose() * (nonlocalEqStrain[0] - localEqStrain[0]) / nonlocalXi[0] +
                                        B.transpose() * (B * mNodalValues.at(Node::NONLOCALEQSTRAIN)));
            break;
        }
        case Node::RELATIVEHUMIDITY:
        {
            const auto& internalGradientRH_B = *static_cast<ConstitutiveVector<TDim>*>(constitutiveOutputMap.at(Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B).get());
            const auto& internalGradientRH_N = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N).get());
            rInternalGradient[dofRow] += mDetJxWeightIPxSection * (mB.at(dofRow).transpose()  * internalGradientRH_B +
                                                                   mN.at(dofRow)->transpose() * internalGradientRH_N);
            break;
        }
        case Node::TEMPERATURE:
        {
            const auto& heatFlux = *static_cast<ConstitutiveVector<TDim>*>(constitutiveOutputMap.at(Constitutive::Output::HEAT_FLUX).get());
            const auto& heatChange = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::HEAT_CHANGE).get());
            rInternalGradient[dofRow] += mDetJxWeightIPxSection * (mB.at(dofRow).transpose() * heatFlux +
                                                                   mN.at(dofRow)->transpose() * heatChange);
            break;
        }
        case Node::WATERVOLUMEFRACTION:
        {
            const auto& internalGradientWV_B = *static_cast<ConstitutiveVector<TDim>*>(constitutiveOutputMap.at(Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B).get());
            const auto& internalGradientWV_N = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N).get());
            rInternalGradient[dofRow] += mDetJxWeightIPxSection * (mB.at(dofRow).transpose()  * internalGradientWV_B +
                                                                   mN.at(dofRow)->transpose() * internalGradientWV_N);
            break;
        }
        case Node::CRACKPHASEFIELD:
        {
            const auto  G       = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::FRACTURE_ENERGY);
            const auto  l       = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER);
            const auto& N       = *(mN.at(Node::CRACKPHASEFIELD));
            const auto& B       = mB.at(Node::CRACKPHASEFIELD);
            const auto& d       = mNodalValues.at(Node::CRACKPHASEFIELD);
            const auto& d_dt    = mNodalValues_dt1.at(Node::CRACKPHASEFIELD);

            const auto& kappa   = (*static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::ELASTIC_ENERGY_DAMAGED_PART).get()))[0];
            const auto  visco   = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::ARTIFICIAL_VISCOSITY);

            //! @todo replace sigma * eps with static data kappa
            rInternalGradient[dofRow] += mDetJxWeightIPxSection * (     ( (G/l + 2.*kappa) * N.transpose() * N
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
void NuTo::ContinuumElement<TDim>::CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0, int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const
{
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
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
            {
                const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
                hessian0 += mDetJxWeightIPxSection *  mB.at(dofRow).transpose() * tangentStressStrain * mB.at(dofRow);
                break;
            }
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN):
            {
                const auto& tangentStressNonlocalEqStrain = *static_cast<ConstitutiveVector<VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN).get());
                hessian0 += mDetJxWeightIPxSection *  mB.at(dofRow).transpose() * tangentStressNonlocalEqStrain * *(mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::DISPLACEMENTS):
            {
                const auto& tangentLocalEqStrainStrain = *static_cast<ConstitutiveVector<VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN).get());
                hessian0 -= mDetJxWeightIPxSection *  mN.at(dofRow)->transpose() * tangentLocalEqStrainStrain.transpose() * (mB.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::NONLOCALEQSTRAIN):
            {
                const auto& N = *(mN.at(dofRow));
                const auto& B = mB.at(dofRow);
                const auto& nonlocalXi = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::NONLOCAL_PARAMETER_XI).get());
                hessian0 += mDetJxWeightIPxSection * (N.transpose() * (1./nonlocalXi[0]) * N + B.transpose() * B);
                break;
            }

            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::TEMPERATURE):
            {
                const auto& tangentHeatFluxTemperatureGradient = *static_cast<ConstitutiveMatrix<TDim, TDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT).get());
                hessian0 += mDetJxWeightIPxSection * mB.at(dofRow).transpose() * tangentHeatFluxTemperatureGradient * mB.at(dofRow);
                break;
            }
                //VHIRTHAMTODO get references to shape functions ---> no double find for the same value
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
            {
                const auto& internalGradientRH_dRH_BB_H0 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0).get());
                const auto& internalGradientRH_dRH_NN_H0 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0).get());
                hessian0 += mDetJxWeightIPxSection * (mB.at(dofRow).transpose()  * internalGradientRH_dRH_BB_H0[0] * mB.at(dofCol) +
                                                      mN.at(dofRow)->transpose() * internalGradientRH_dRH_NN_H0 * (*mN.at(dofCol)));
                break;
            }
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
            {
                const auto& internalGradientRH_dWV_BN_H0 = *static_cast<ConstitutiveVector<TDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0).get());
                const auto& internalGradientRH_dWV_NN_H0 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0).get());
                hessian0 += mDetJxWeightIPxSection * (mB.at(dofRow).transpose()  * internalGradientRH_dWV_BN_H0 * (*mN.at(dofCol)) +
                                                      mN.at(dofRow)->transpose() * internalGradientRH_dWV_NN_H0 * (*mN.at(dofCol)));
                break;
            }
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
            {
                const auto& internalGradientWV_dRH_NN_H0 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0).get());
                hessian0 += mDetJxWeightIPxSection * mN.at(dofRow)->transpose()  * internalGradientWV_dRH_NN_H0 * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
            {
                const auto& internalGradientWV_dWV_BB_H0 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0).get());
                const auto& internalGradientWV_dWV_BN_H0 = *static_cast<ConstitutiveVector<TDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0).get());
                const auto& internalGradientWV_dWV_NN_H0 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0).get());
                hessian0 += mDetJxWeightIPxSection * (mB.at(dofRow).transpose()  * internalGradientWV_dWV_BB_H0[0] * mB.at(dofCol) +
                                                      mB.at(dofRow).transpose()  * internalGradientWV_dWV_BN_H0 * (*mN.at(dofCol)) +
                                                      mN.at(dofRow)->transpose() * internalGradientWV_dWV_NN_H0 * (*mN.at(dofCol)));
                break;
            }
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::RELATIVEHUMIDITY):
            {
                const auto& engineeringStress_dRH = *static_cast<ConstitutiveVector<VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY).get());
                hessian0 += mDetJxWeightIPxSection * mB.at(dofRow).transpose()  * engineeringStress_dRH * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE):
            {
                const auto& dStressDTemperature = *static_cast<ConstitutiveVector<VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_ENGINEERING_STRESS_D_TEMPERATURE).get());
                hessian0 += mDetJxWeightIPxSection * mB.at(dofRow).transpose()  * dStressDTemperature * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::WATERVOLUMEFRACTION):
            {
                const auto& engineeringStress_dWV = *static_cast<ConstitutiveVector<VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION).get());
                hessian0 += mDetJxWeightIPxSection * mB.at(dofRow).transpose()  * engineeringStress_dWV * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::CRACKPHASEFIELD):
            {
                const auto& Bu = mB.at(Node::eDof::DISPLACEMENTS);
                const auto& Nd = *(mN.at(Node::eDof::CRACKPHASEFIELD));
                const auto& dStressDPhaseField = *static_cast<ConstitutiveVector<VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_ENGINEERING_STRESS_D_PHASE_FIELD).get());
                hessian0 += mDetJxWeightIPxSection * Bu.transpose() * dStressDPhaseField * Nd;
                break;
            }

            case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::DISPLACEMENTS):
            {
                const auto& Nd = *(mN.at(Node::eDof::CRACKPHASEFIELD));
                const double d = (Nd * mNodalValues.at(Node::eDof::CRACKPHASEFIELD)).at(0,0);

                const auto& Bu = mB.at(Node::eDof::DISPLACEMENTS);
                const auto& tangentElasticEnergyStrain = *static_cast<ConstitutiveVector<VoigtDim>*>(constitutiveOutputMap.at(Constitutive::Output::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN).get());

                hessian0 += mDetJxWeightIPxSection * 2 * (d -1.) * Nd.transpose() * tangentElasticEnergyStrain.transpose() * Bu;
                break;
            }

            case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::CRACKPHASEFIELD):
            {
                const auto  G = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::FRACTURE_ENERGY);
                const auto  l = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER);

                const auto& N = *(mN.at(Node::eDof::CRACKPHASEFIELD));
                const auto& B = mB.at(Node::eDof::CRACKPHASEFIELD);

                const auto& kappa   = (*static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::ELASTIC_ENERGY_DAMAGED_PART).get()))[0];

                hessian0 += mDetJxWeightIPxSection * ( (G/l + 2*kappa) * N.transpose() * N +  G*l*B.transpose() * B );
                break;
            }
            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::DISPLACEMENTS):
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for (" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputHessian1(BlockFullMatrix<double>& rHessian1, int rTheIP,
        const ConstitutiveOutputMap& constitutiveOutputMap) const
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
            {
                const auto& tangentHeatTemperature = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_HEAT_D_TEMPERATURE).get());
                hessian1 += mDetJxWeightIPxSection * mN.at(dofRow)->transpose()
                          * tangentHeatTemperature * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
            {
                const auto& internalGradientRH_dRH_NN_H1 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1).get());
                hessian1 += mDetJxWeightIPxSection * mN.at(dofRow)->transpose() * internalGradientRH_dRH_NN_H1 * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
            {
                const auto& internalGradientRH_dWV_NN_H1 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1).get());
                hessian1 += mDetJxWeightIPxSection * mN.at(dofRow)->transpose() * internalGradientRH_dWV_NN_H1 * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
                break;

            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
            {
                const auto& internalGradientWV_dWV_NN_H1 = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1).get());
                hessian1 += mDetJxWeightIPxSection * mN.at(dofRow)->transpose() * internalGradientWV_dWV_NN_H1 * (*mN.at(dofCol));
                break;
            }
            case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::CRACKPHASEFIELD):
            {
                const auto& N = *(mN.at(Node::eDof::CRACKPHASEFIELD));
                const auto  visco   = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::ARTIFICIAL_VISCOSITY);
                hessian1 += visco * mDetJxWeightIPxSection * N.transpose() * N;
            }
                break;
            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::TEMPERATURE, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE):
                break;
            default:
                throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Element output HESSIAN_1_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumElement<TDim>::CalculateElementOutputHessian2(BlockFullMatrix<double>& rHessian2, int rTheIP) const
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
                const auto& N = *(mN.at(dofRow));
                double rho = GetConstitutiveLaw(rTheIP)->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY);
                hessian2 += rho * N.transpose() * N * mDetJxWeightIPxSection;

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
void NuTo::ContinuumElement<TDim>::CalculateElementOutputIpData(ElementOutputIpData& rIpData,
        int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const
{
    for (auto& it : rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    {                                       // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::DAMAGE:
            it.second.col(rTheIP) = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::DAMAGE).get());
            break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
            it.second.col(rTheIP) = *static_cast<EngineeringStrain<3>*>(constitutiveOutputMap.at(Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_VISUALIZE).get());
            break;
        case NuTo::IpData::ENGINEERING_STRAIN:
            it.second.col(rTheIP) = *static_cast<EngineeringStrain<3>*>(constitutiveOutputMap.at(Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE).get());
            break;
        case NuTo::IpData::ENGINEERING_STRESS:
            it.second.col(rTheIP) = *static_cast<EngineeringStress<3>*>(constitutiveOutputMap.at(Constitutive::Output::ENGINEERING_STRESS_VISUALIZE).get());
            break;
        case NuTo::IpData::EXTRAPOLATION_ERROR:
            it.second.col(rTheIP) = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::EXTRAPOLATION_ERROR).get());
            break;
        case NuTo::IpData::LOCAL_EQ_STRAIN:
            it.second.col(rTheIP) = *static_cast<ConstitutiveScalar*>(constitutiveOutputMap.at(Constitutive::Output::LOCAL_EQ_STRAIN).get());
            break;
        case NuTo::IpData::SHRINKAGE_STRAIN:
            it.second.col(rTheIP) = *static_cast<EngineeringStrain<3>*>(constitutiveOutputMap.at(Constitutive::Output::SHRINKAGE_STRAIN_VISUALIZE).get());
            break;
        case NuTo::IpData::THERMAL_STRAIN:
            it.second.col(rTheIP) = *static_cast<EngineeringStrain<3>*>(constitutiveOutputMap.at(Constitutive::Output::THERMAL_STRAIN).get());
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
void NuTo::ContinuumElement<TDim>::CalculateNMatrixBMatrixDetJacobian(int rTheIP)
{
    // calculate Jacobian
    const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(rTheIP);

    Eigen::Matrix<double, TDim, TDim> jacobian = CalculateJacobian(derivativeShapeFunctionsGeometryNatural, mNodalValues.at(Node::COORDINATES));
    mDetJacobian = jacobian.determinant();
    if (mDetJacobian == 0)
    {
        for (auto* node : mNodes)
        {
            std::cout << "ElementID " << mStructure->ElementGetId(this) << std::endl;
            std::cout << "NodeID: " << mStructure->NodeGetId(node) << std::endl;
            std::cout << node->Get(Node::COORDINATES) << std::endl << std::endl;
        }
        std::cout << mNodalValues.at(Node::COORDINATES) << std::endl;
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Determinant of the Jacobian is zero, no inversion possible.");
    }

    Eigen::Matrix<double, TDim, TDim> invJacobian = jacobian.inverse();

    // calculate shape functions and their derivatives
    for (auto dof : mInterpolationType->GetDofs())
    {
        if (dof == Node::COORDINATES)
            continue;
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        mN[dof] = &interpolationType.GetMatrixN(rTheIP);

        mB[dof] = CalculateMatrixB(dof, interpolationType.GetDerivativeShapeFunctionsNatural(rTheIP), invJacobian);
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
double NuTo::ContinuumElement<1>::CalculateDetJxWeightIPxSection(int rTheIP) const
{
    Eigen::MatrixXd matrixN = mInterpolationType->Get(Node::COORDINATES).GetMatrixN(rTheIP);
    Eigen::VectorXd globalIPCoordinate = matrixN * ExtractNodeValues(0, Node::COORDINATES);

    return mDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetArea() * mSection->AsSectionTruss()->GetAreaFactor(globalIPCoordinate.at(0, 0));
}

template<>
double NuTo::ContinuumElement<2>::CalculateDetJxWeightIPxSection(int rTheIP) const
{
    return mDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) * mSection->GetThickness();
}

template<>
double NuTo::ContinuumElement<3>::CalculateDetJxWeightIPxSection(int rTheIP) const
{
    return mDetJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP);
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
