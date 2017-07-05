#include "mechanics/constitutive/staticData/IPAdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"

using namespace NuTo::Constitutive;


IPAdditiveInputExplicit::IPAdditiveInputExplicit(AdditiveInputExplicit& rLaw)
    : mLaw(rLaw)
{
    mMainLawIP = mLaw.mMainLaw->CreateIPLaw();
    for (auto& sublaw : mLaw.mSublaws)
    {
        mSublawIPs.push_back(sublaw->CreateIPLaw().release());
    }
}


IPAdditiveInputExplicit::IPAdditiveInputExplicit(const IPAdditiveInputExplicit& rOther)
    : mLaw(rOther.mLaw)
    , mSublawIPs(rOther.mSublawIPs)
{
    this->mMainLawIP = rOther.mMainLawIP->Clone();
}


std::unique_ptr<IPConstitutiveLawBase> IPAdditiveInputExplicit::Clone() const
{
    return std::make_unique<IPAdditiveInputExplicit>(*this);
}


NuTo::ConstitutiveBase& IPAdditiveInputExplicit::GetConstitutiveLaw() const
{
    return mLaw;
}


void IPAdditiveInputExplicit::AllocateAdditional(int rNum)
{
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.AllocateAdditional(rNum);
    }
}


void IPAdditiveInputExplicit::ShiftToPast()
{
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.ShiftToPast();
    }
}


void IPAdditiveInputExplicit::ShiftToFuture()
{
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.ShiftToFuture();
    }
}


void IPAdditiveInputExplicit::NuToSerializeSave(SerializeStreamOut& rStream)
{
    IPConstitutiveLawBase::NuToSerializeSave(rStream);
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.NuToSerializeSave(rStream);
    }
}


void IPAdditiveInputExplicit::NuToSerializeLoad(SerializeStreamIn& rStream)
{
    IPConstitutiveLawBase::NuToSerializeLoad(rStream);
    for (auto& sublaw : mSublawIPs)
    {
        sublaw.NuToSerializeLoad(rStream);
    }
}


template <int TDim>
void IPAdditiveInputExplicit::AdditiveInputExplicitEvaluate(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                            const NuTo::ConstitutiveOutputMap& rConstitutiveOutput)
{
    // Copy inputs for main law, because they might be modified by the sublaws and these modifications will be
    // passed above the borders of this law.
    NuTo::ConstitutiveInputMap mainLawInputMap = rConstitutiveInput;

    // evaluate sublaws
    std::vector<NuTo::ConstitutiveOutputMap> sublawOutputMapVec;

    for (size_t i = 0; i < mSublawIPs.size(); ++i)
    {
        // Get the sublaw specific output map depending on the main laws inputs and the global outputs
        sublawOutputMapVec.emplace_back(mLaw.GetSublawOutputMap<TDim>(rConstitutiveInput, rConstitutiveOutput, i));

        mSublawIPs[i].Evaluate<TDim>(rConstitutiveInput, sublawOutputMapVec[i]);

        // Apply outputs to the main laws input and the global outputs
        ApplySublawOutputs<TDim>(mainLawInputMap, rConstitutiveOutput, sublawOutputMapVec[i]);
    }

    mMainLawIP->Evaluate<TDim>(mainLawInputMap, rConstitutiveOutput);

    // calculate derivatives that depend on outputs from the main law and the sublaws
    CalculateDerivatives<TDim>(rConstitutiveOutput, sublawOutputMapVec);
}


template <int TDim>
void IPAdditiveInputExplicit::CalculateDerivatives(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                   std::vector<ConstitutiveOutputMap>& rSublawOutputVec)
{
    constexpr const int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
        {
            for (size_t i = 0; i < mSublawIPs.size(); ++i)
            {
                if (mLaw.mInputsToModify[i] == Constitutive::eInput::ENGINEERING_STRAIN)
                {
                    // Get the corresponding sublaw output
                    Constitutive::eOutput derivativeOutputEnum =
                            mLaw.GetDerivativeEnumSublaw(Constitutive::eOutput::ENGINEERING_STRAIN, itOutput.first);

                    const ConstitutiveOutputMap::iterator& sublawOutput =
                            rSublawOutputVec[i].find(derivativeOutputEnum);

                    if (sublawOutput == rSublawOutputVec[i].end())
                        // if current sublaw, does not provide the needed output, continue with next law
                        continue;
                    else
                    {
                        assert(rConstitutiveOutput.count(
                                Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));
                        assert(itOutput.second->GetIsCalculated() == false &&
                               "Currently, it is not supported that multiple sublaws write to the same derivative.");
                        if (sublawOutput->second->GetIsCalculated() == false)
                            throw MechanicsException(
                                    __PRETTY_FUNCTION__,
                                    "The value " + Constitutive::OutputToString(sublawOutput->first) +
                                            ", which is necessary to determine " +
                                            Constitutive::OutputToString(itOutput.first) +
                                            " was requested from a sublaw but has not been calculated!");
                        const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(
                                rConstitutiveOutput.at(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)
                                        .get());
                        const auto& sublawDerivative =
                                *static_cast<ConstitutiveVector<VoigtDim>*>(sublawOutput->second.get());

                        static_cast<EngineeringStress<TDim>*>(itOutput.second.get())->AsVector() =
                                tangentStressStrain * sublawDerivative;
                    }
                }
                else
                    // if law does not modify strains, continue with next sublaw
                    continue;
            }
            break;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}


template <int TDim>
void IPAdditiveInputExplicit::ApplySublawOutputs(const ConstitutiveInputMap& rMainLawInput,
                                                 const ConstitutiveOutputMap& rConstitutiveOutput,
                                                 const ConstitutiveOutputMap& rSublawOutput)
{
    for (const auto& it : rSublawOutput)
    {
        if (it.first == Constitutive::eOutput::ENGINEERING_STRAIN)
        {
            // Modify input strain for main constitutive law
            assert(rSublawOutput.at(Constitutive::eOutput::ENGINEERING_STRAIN) != nullptr);
            assert(rMainLawInput.at(Constitutive::eInput::ENGINEERING_STRAIN) != nullptr);
            *static_cast<EngineeringStrain<TDim>*>(rMainLawInput.at(Constitutive::eInput::ENGINEERING_STRAIN).get()) -=
                    *static_cast<EngineeringStrain<TDim>*>(
                            rSublawOutput.at(Constitutive::eOutput::ENGINEERING_STRAIN).get());
        }
        else if (rConstitutiveOutput.count(it.first) and rConstitutiveOutput.at(it.first) != nullptr)
        {
            // Copy sublaw output data direct into global output
            *(rConstitutiveOutput.at(it.first)) = *(it.second);
        }
    }
}


template void
IPAdditiveInputExplicit::AdditiveInputExplicitEvaluate<1>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                          const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
template void
IPAdditiveInputExplicit::AdditiveInputExplicitEvaluate<2>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                          const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
template void
IPAdditiveInputExplicit::AdditiveInputExplicitEvaluate<3>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                          const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
