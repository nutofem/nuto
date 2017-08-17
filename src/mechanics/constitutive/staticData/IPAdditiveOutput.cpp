#include "mechanics/constitutive/staticData/IPAdditiveOutput.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

using namespace NuTo::Constitutive;


IPAdditiveOutput::IPAdditiveOutput(AdditiveOutput& rLaw)
    : mLaw(rLaw)
{
    for (auto& sublaw : mLaw.mSublaws)
    {
        mSublawIPs.push_back(sublaw->CreateIPLaw().release());
    }
}


IPAdditiveOutput::IPAdditiveOutput(const IPAdditiveOutput& rOther)
    : IPConstitutiveLawBase(rOther)
    , mLaw(rOther.mLaw)
    , mSublawIPs(rOther.mSublawIPs)
{
}


std::unique_ptr<IPConstitutiveLawBase> IPAdditiveOutput::Clone() const
{
    return std::make_unique<IPAdditiveOutput>(*this);
}


NuTo::ConstitutiveBase& IPAdditiveOutput::GetConstitutiveLaw() const
{
    return mLaw;
}


void IPAdditiveOutput::AllocateAdditional(int rNum)
{
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.AllocateAdditional(rNum);
    }
}


void IPAdditiveOutput::ShiftToPast()
{
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.ShiftToPast();
    }
}


void IPAdditiveOutput::ShiftToFuture()
{
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.ShiftToFuture();
    }
}


void IPAdditiveOutput::NuToSerializeSave(SerializeStreamOut& rStream)
{
    IPConstitutiveLawBase::NuToSerializeSave(rStream);
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.NuToSerializeSave(rStream);
    }
}


void IPAdditiveOutput::NuToSerializeLoad(SerializeStreamIn& rStream)
{
    IPConstitutiveLawBase::NuToSerializeLoad(rStream);
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.NuToSerializeLoad(rStream);
    }
}


IPConstitutiveLawBase* IPAdditiveOutput::GetSublawIP(ConstitutiveBase* rCLPtr)
{
    for (auto& sublawIP : mSublawIPs)
    {
        ConstitutiveBase& sublaw = sublawIP.GetConstitutiveLaw();
        if ((&sublaw) == rCLPtr)
            return (&sublawIP);
        switch (sublaw.GetType())
        {
        case Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT:
        case Constitutive::eConstitutiveType::ADDITIVE_INPUT_IMPLICIT:
        case Constitutive::eConstitutiveType::ADDITIVE_OUTPUT:
        {
            IPConstitutiveLawBase* additiveSublawOutput = sublawIP.GetSublawIP(rCLPtr);
            if (additiveSublawOutput != nullptr)
                return additiveSublawOutput;
            break;
        }
        default:
            break;
        }
    }
    return nullptr;
}


template <int TDim>
void IPAdditiveOutput::AdditiveOutputEvaluate(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                              const NuTo::ConstitutiveOutputMap& rConstitutiveOutput)
{
    using namespace Constitutive;
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    for (auto& output : rConstitutiveOutput)
    {
        if (output.second != nullptr)
            output.second->SetZero();
    }

    for (auto& sublaw : mSublawIPs)
    {
        NuTo::ConstitutiveOutputMap singleOutput;
        for (auto& output : rConstitutiveOutput)
        {
            singleOutput[output.first] = ConstitutiveIOBase::makeConstitutiveIO<TDim>(output.first);
        }

        sublaw.Evaluate<TDim>(rConstitutiveInput, singleOutput);

        for (const auto& output : singleOutput)
        {
            if (output.second != nullptr and output.second->GetIsCalculated())
            {
                switch (output.first)
                {
                case eOutput::LOCAL_EQ_STRAIN:
                case eOutput::DAMAGE:
                case eOutput::EXTRAPOLATION_ERROR:
                case eOutput::HEAT_CHANGE:
                case eOutput::D_HEAT_D_TEMPERATURE:
                case eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
                case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
                case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
                case eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
                case eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
                case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
                case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
                case eOutput::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
                case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
                case eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
                case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
                case eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
                case eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
                case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
                case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
                case eOutput::NONLOCAL_RADIUS:
                    *static_cast<ConstitutiveScalar*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveScalar*>(singleOutput.at(output.first).get());
                    rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                    break;
                case eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
                case eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
                case eOutput::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
                case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
                case eOutput::HEAT_FLUX:
                    *static_cast<ConstitutiveVector<TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveVector<TDim>*>(singleOutput.at(output.first).get());
                    rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                    break;
                case eOutput::ENGINEERING_STRESS:
                case eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
                case eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
                case eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
                case eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
                case eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN:
                    *static_cast<ConstitutiveVector<VoigtDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveVector<VoigtDim>*>(singleOutput.at(output.first).get());
                    rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                    break;
                case eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
                case eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1:
                    *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(singleOutput.at(output.first).get());
                    rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                    break;
                case eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
                case eOutput::ENGINEERING_STRAIN_VISUALIZE:
                case eOutput::SHRINKAGE_STRAIN_VISUALIZE:
                case eOutput::THERMAL_STRAIN:
                    *static_cast<EngineeringStrain<TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<EngineeringStrain<TDim>*>(singleOutput.at(output.first).get());
                    rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                    break;
                case eOutput::ENGINEERING_STRESS_VISUALIZE:
                    *static_cast<EngineeringStress<TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<EngineeringStress<TDim>*>(singleOutput.at(output.first).get());
                    rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                    break;
                case eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
                    *static_cast<ConstitutiveMatrix<TDim, TDim>*>(rConstitutiveOutput.at(output.first).get()) +=
                            *static_cast<ConstitutiveMatrix<TDim, TDim>*>(singleOutput.at(output.first).get());
                    rConstitutiveOutput.at(output.first)->SetIsCalculated(true);
                    break;

                default:
                    throw Exception(__PRETTY_FUNCTION__, "Output is not implemented or can't be handled.");
                } // switch outputs
            } // if not nullptr and IsCalculated
        } // for each output
    } // for each sublaw
}


template void IPAdditiveOutput::AdditiveOutputEvaluate<1>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                          const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
template void IPAdditiveOutput::AdditiveOutputEvaluate<2>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                          const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
template void IPAdditiveOutput::AdditiveOutputEvaluate<3>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                          const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
