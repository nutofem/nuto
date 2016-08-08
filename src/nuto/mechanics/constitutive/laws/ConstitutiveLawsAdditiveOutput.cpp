#include "ConstitutiveLawsAdditiveOutput.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"

//VHIRTHAMTODO --- maybe store possible dof combinations one time in a set when constitutive law is added instead of looping over each law every time
bool NuTo::ConstitutiveLawsAdditiveOutput::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    if(mComputableDofCombinations[rTimeDerivative].find(std::pair<Node::eDof,Node::eDof>(rDofRow,rDofCol)) != mComputableDofCombinations[rTimeDerivative].end())
        return true;
    return false;
}

template <int TDim>
NuTo::Error::eError NuTo::ConstitutiveLawsAdditiveOutput::Evaluate(NuTo::ElementBase* rElement, int rIp,
        const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    using namespace Constitutive::Output;
    Error::eError error = Error::SUCCESSFUL;
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    for (auto& output : rConstitutiveOutput)
    {
        if (output.second != nullptr) output.second->SetZero();
    }

    try
    {
        for(unsigned int i = 0; i < mConstitutiveLaws.size(); ++i)
        {
            NuTo::ConstitutiveOutputMap singleOutput;
            for (auto& output : rConstitutiveOutput)
            {
                singleOutput[output.first] = ConstitutiveIOBase::makeConstitutiveIO<TDim>(output.first);
            }
            if (TDim == 1)
            {
                error = mConstitutiveLaws[i]->Evaluate1D(rElement, rIp, rConstitutiveInput, singleOutput);
            }
            if (TDim == 2)
            {
                error = mConstitutiveLaws[i]->Evaluate2D(rElement, rIp, rConstitutiveInput, singleOutput);
            }
            if (TDim == 3)
            {
                error = mConstitutiveLaws[i]->Evaluate3D(rElement, rIp, rConstitutiveInput, singleOutput);
            }
            for (const auto& output : singleOutput)
            {
                if (output.second != nullptr and output.second->GetIsCalculated())
                {
                    switch (output.first)
                    {
                    case LOCAL_EQ_STRAIN:
                    case NONLOCAL_PARAMETER_XI:
                    case DAMAGE:
                    case EXTRAPOLATION_ERROR:
                    case HEAT_CHANGE:
                    case D_HEAT_D_TEMPERATURE:
                    case INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
                    case D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
                    case D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
                    case D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
                    case INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
                    case D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
                    case D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
                    case D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
                    case D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
                    case D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
                    case D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
                    case INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
                    case INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
                    case D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
                    case D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
                        *static_cast<ConstitutiveScalar*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveScalar*>(singleOutput.at(output.first).get());
                        rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                        break;
                    case INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
                    case INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
                    case D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
                    case D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
                    case HEAT_FLUX:
                        *static_cast<ConstitutiveVector<TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveVector<TDim>*>(singleOutput.at(output.first).get());
                        rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                        break;
                    case ENGINEERING_STRESS:
                    case D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
                    case D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
                    case D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
                    case D_LOCAL_EQ_STRAIN_XI_D_STRAIN:
                    case D_ENGINEERING_STRESS_D_TEMPERATURE:
                    case D_LOCAL_EQ_STRAIN_D_STRAIN:
                        *static_cast<ConstitutiveVector<VoigtDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveVector<VoigtDim>*>(singleOutput.at(output.first).get());
                        rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                        break;
                    case D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
                        *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(singleOutput.at(output.first).get());
                        rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                        break;
                    case ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
                    case ENGINEERING_STRAIN_VISUALIZE:
                    case SHRINKAGE_STRAIN_VISUALIZE:
                    case THERMAL_STRAIN:
                        *static_cast<EngineeringStrain<TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<EngineeringStrain<TDim>*>(singleOutput.at(output.first).get());
                        rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                        break;
                    case ENGINEERING_STRESS_VISUALIZE:
                        *static_cast<EngineeringStress<TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<EngineeringStress<TDim>*>(singleOutput.at(output.first).get());
                        rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                        break;
                    case D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
                        *static_cast<ConstitutiveMatrix<TDim, TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveMatrix<TDim, TDim>*>(singleOutput.at(output.first).get());
                        rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                        break;
                    } // switch outputs
                } // if not nullptr and IsCalculated
            } // for each output
        } // for each sublaw
    } //try
    catch(Exception e)
    {
        e.AddMessage(__PRETTY_FUNCTION__,"Exception while evaluating constitutive law attached to an additive output.");
        throw;
    }
    return error;

}


NuTo::Error::eError NuTo::ConstitutiveLawsAdditiveOutput::Evaluate1D(NuTo::ElementBase *rElement, int rIp,
        const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    return NuTo::ConstitutiveLawsAdditiveOutput::Evaluate<1>(rElement, rIp, rConstitutiveInput, rConstitutiveOutput);
}


NuTo::Error::eError NuTo::ConstitutiveLawsAdditiveOutput::Evaluate2D(NuTo::ElementBase *rElement, int rIp,
        const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    return NuTo::ConstitutiveLawsAdditiveOutput::Evaluate<2>(rElement, rIp, rConstitutiveInput, rConstitutiveOutput);
}


NuTo::Error::eError NuTo::ConstitutiveLawsAdditiveOutput::Evaluate3D(NuTo::ElementBase *rElement, int rIp,
        const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    return NuTo::ConstitutiveLawsAdditiveOutput::Evaluate<3>(rElement, rIp, rConstitutiveInput, rConstitutiveOutput);
}


NuTo::ConstitutiveInputMap NuTo::ConstitutiveLawsAdditiveOutput::GetConstitutiveInputs(
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput, const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
    {
        ConstitutiveInputMap singleLawInputMap = mConstitutiveLaws[i]->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                             rInterpolationType);
        constitutiveInputMap.Merge(singleLawInputMap);
    }

    return constitutiveInputMap;
}


void NuTo::ConstitutiveLawsAdditiveOutput::AddCalculableDofCombinations(NuTo::ConstitutiveBase *rConstitutiveLaw)
{
    std::set<Node::eDof> allDofs = Node::GetDofSet();
    for (unsigned int i=0; i<mComputableDofCombinations.size(); ++i)
        for (auto itRow : allDofs)
            for (auto itCol : allDofs)
            {
                if (rConstitutiveLaw->CheckDofCombinationComputable(itRow,itCol,i))
                        mComputableDofCombinations[i].emplace(itRow,itCol);
            }
}
