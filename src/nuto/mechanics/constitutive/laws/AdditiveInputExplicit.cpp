#include "AdditiveInputExplicit.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/staticData/Composite.h"

template <int TDim>
NuTo::Constitutive::StaticData::Component* NuTo::AdditiveInputExplicit::AllocateStaticData(const NuTo::ElementBase *rElement) const
{
    auto composite =
        static_cast<Constitutive::StaticData::Composite*>(AdditiveBase::AllocateStaticData<TDim>(rElement));
    Constitutive::StaticData::Component* mainComponent;
    switch (TDim)
    {
    case 1:
        mainComponent = mMainLaw->AllocateStaticData1D(rElement);
        break;
    case 2:
        mainComponent = mMainLaw->AllocateStaticData2D(rElement);
        break;
    case 3:
        mainComponent = mMainLaw->AllocateStaticData3D(rElement);
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Invalid dimension.");
    }
    composite->AddComponent(mainComponent);
    return composite;
}

void NuTo::AdditiveInputExplicit::AddConstitutiveLaw(NuTo::ConstitutiveBase& rConstitutiveLaw,
        Constitutive::eInput rModiesInput)
{
    if(rModiesInput == Constitutive::eInput::NONE)
    {
        if (mMainLaw != nullptr)
            throw MechanicsException(__PRETTY_FUNCTION__,
                "There can be only one! --- This additive input law only accepts one law which calculates the "
                "output. All other laws are only allowed to modify the input to this law. Specify the modifying "
                "laws by providing the enum of the modified input as second function parameter.");
        mMainLaw = &rConstitutiveLaw;
        AddCalculableDofCombinations(rConstitutiveLaw);
    }
    else
    {
        mInputsToModify.push_back(rModiesInput);
        AdditiveBase::AddConstitutiveLaw(rConstitutiveLaw, rModiesInput); 
    }
}


NuTo::ConstitutiveInputMap NuTo::AdditiveInputExplicit::GetConstitutiveInputs(
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
        const NuTo::InterpolationType &rInterpolationType) const
{
    // Get Inputs for output returning constitutive law
    ConstitutiveInputMap mainLawConstitutiveInputMap(mMainLaw->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                     rInterpolationType));
    // Get Inputs for input modifying constitutive laws
    ConstitutiveInputMap sublawsConstitutiveInputMap;
    for (unsigned int i = 0; i < mSublaws.size(); ++i)
    {

        // Get the necessary inputs for the sublaws ---> the modifications to the main laws inputs are returned as
        // output from the sublaws. Therefore one first needs the sublaws output list, depending on the main laws
        // inputs and the global outputs.

        // INFO regarding the template function call: The dimension can be chosen freely, because the created object is
        // never used. Should be 1 because of the lowest construction costs.  Alternative would be to copy the called
        // function and replace all objects with nullptr. But then you have to maintain 2 nearly identical functions.
        ConstitutiveOutputMap sublawOutputMap = GetSublawOutputMap<1>(mainLawConstitutiveInputMap,
                rConstitutiveOutput, i);

        // Don't merge the sublaw inputs directly into the main laws input map!  ---> When more than one sublaw is
        // attached the inputs of the first might effect the following laws outputs - have a look at the line above!
        sublawsConstitutiveInputMap.Merge(mSublaws[i]->GetConstitutiveInputs(sublawOutputMap, rInterpolationType));
    }
    return mainLawConstitutiveInputMap.Merge(sublawsConstitutiveInputMap);
}

NuTo::Constitutive::eConstitutiveType NuTo::AdditiveInputExplicit::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT;
}


template <int TDim>
void NuTo::AdditiveInputExplicit::ApplySublawOutputs(const ConstitutiveInputMap& rMainLawInput,
        const ConstitutiveOutputMap& rConstitutiveOutput, const ConstitutiveOutputMap& rSublawOutput)
{
    for (const auto& it : rSublawOutput)
    {
        if (it.first == Constitutive::eOutput::ENGINEERING_STRAIN)
        {
            // Modify input strain for main constitutive law
            assert(rSublawOutput.at(Constitutive::eOutput::ENGINEERING_STRAIN)!=nullptr);
            assert(rMainLawInput.at(Constitutive::eInput::ENGINEERING_STRAIN)!=nullptr);
            *static_cast<EngineeringStrain<TDim>*>(rMainLawInput.at(Constitutive::eInput::ENGINEERING_STRAIN).get()) -= *static_cast<EngineeringStrain<TDim>*>(rSublawOutput.at(Constitutive::eOutput::ENGINEERING_STRAIN).get());
        }
        else if (rConstitutiveOutput.count(it.first) and rConstitutiveOutput.at(it.first) != nullptr)
        {
            //Copy sublaw output data direct into global output
            *(rConstitutiveOutput.at(it.first)) = *(it.second);
        }
    }
}


template <int TDim>
void NuTo::AdditiveInputExplicit::CalculateDerivatives(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                       std::vector<ConstitutiveOutputMap>& rSublawOutputVec)
{
    constexpr const int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
        {
            for (unsigned int i = 0; i < mSublaws.size(); ++i)
            {
                if(mInputsToModify[i] == Constitutive::eInput::ENGINEERING_STRAIN)
                {
                    // Get the corresponding sublaw output
                    Constitutive::eOutput derivativeOutputEnum =
                        GetDerivativeEnumSublaw(Constitutive::eOutput::ENGINEERING_STRAIN, itOutput.first);

                    const ConstitutiveOutputMap::iterator& sublawOutput = rSublawOutputVec[i].find(derivativeOutputEnum);

                    if(sublawOutput == rSublawOutputVec[i].end())
                        // if current sublaw, does not provide the needed output, continue with next law
                        continue;
                    else
                    {
                        assert(rConstitutiveOutput.count(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));
                        assert(itOutput.second->GetIsCalculated() == false &&
                                "Currently it is not supported that multiple sublaws write to the same derivative");
                        if(sublawOutput->second->GetIsCalculated() == false)
                            throw MechanicsException(__PRETTY_FUNCTION__, "The value "
                                    + Constitutive::OutputToString(sublawOutput->first) + 
                                    ", which is necessary to determine " + Constitutive::OutputToString(itOutput.first)
                                    + " was requested  from a sublaw but has not been calculated!" );
                        const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<VoigtDim, VoigtDim>*>(rConstitutiveOutput.at(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
                        const auto& sublawDerivative = *static_cast<ConstitutiveVector<VoigtDim>*>(sublawOutput->second.get());

                        (static_cast<EngineeringStress<TDim>*>(itOutput.second.get()))->AsVector() = tangentStressStrain * sublawDerivative;
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
NuTo::eError NuTo::AdditiveInputExplicit::EvaluateAdditiveInputExplicit(
        const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
        NuTo::Constitutive::StaticData::Component* staticData)
{
    eError error = eError::SUCCESSFUL;
    // Copy inputs for main law, because they might be modified by the sublaws and these modifications will be passed above the borders of this law.
    NuTo::ConstitutiveInputMap mainLawInputMap = rConstitutiveInput;

    // evaluate sublaws
    std::vector<NuTo::ConstitutiveOutputMap> sublawOutputMapVec;

    auto& localStaticData = *dynamic_cast<Constitutive::StaticData::Composite*>(staticData);
    for (unsigned int i = 0; i < mSublaws.size(); ++i)
    {
        // Get the sublaw specific output map depending on the main laws inputs and the global outputs
        sublawOutputMapVec.emplace_back(GetSublawOutputMap<TDim>(rConstitutiveInput, rConstitutiveOutput, i));

        mSublaws[i]->Evaluate<TDim>(rConstitutiveInput, sublawOutputMapVec[i], &localStaticData.GetComponent(i+1));

        if(error != eError::SUCCESSFUL)
            throw MechanicsException(__PRETTY_FUNCTION__,
                    "Attached constitutive law returned an error code. Can't handle this");

        // Apply outputs to the main laws input and the global outputs
        ApplySublawOutputs<TDim>(mainLawInputMap, rConstitutiveOutput, sublawOutputMapVec[i]);
    }

    error = mMainLaw->Evaluate<TDim>(mainLawInputMap, rConstitutiveOutput, &localStaticData.GetComponent(0));

    // calculate derivatives that depend on outputs from the main law and the sublaws
    CalculateDerivatives<TDim>(rConstitutiveOutput, sublawOutputMapVec);

    return error;
}


NuTo::Constitutive::eOutput NuTo::AdditiveInputExplicit::GetDerivativeEnumSublaw(NuTo::Constitutive::eOutput rParameter,
                                                                                         NuTo::Constitutive::eOutput rMainDerivative) const
{
    switch (rParameter)
    {
    case Constitutive::eOutput::ENGINEERING_STRAIN:
        switch(rMainDerivative)
        {
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
            return Constitutive::eOutput::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY;

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
            return Constitutive::eOutput::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION;

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
            return Constitutive::eOutput::D_STRAIN_D_TEMPERATURE;

        default:
            throw Exception(__PRETTY_FUNCTION__,
                            std::string("No partial derivative defined for parameter " + Constitutive::OutputToString(rParameter))
                            + " and global derivative " + Constitutive::OutputToString(rMainDerivative));
        }
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        std::string("No partial derivatives defined for parameter " + Constitutive::OutputToString(rParameter)));
    }
}


template <int TDim>
NuTo::ConstitutiveOutputMap NuTo::AdditiveInputExplicit::GetSublawOutputMap(
        const NuTo::ConstitutiveInputMap& rMainLawInputMap, const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
        unsigned int rSublawIndex) const
{
    ConstitutiveOutputMap modifiedOutputMap;

    // Add sublaw outputs depending on main law outputs
    for (const auto& itMainLawOutput : rMainLawOutputMap)
    {
        switch (itMainLawOutput.first)
        {
        case Constitutive::eOutput::SHRINKAGE_STRAIN_VISUALIZE:
        case Constitutive::eOutput::THERMAL_STRAIN:
            modifiedOutputMap.emplace(itMainLawOutput.first,itMainLawOutput.second->clone());
            break;
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
            if(mInputsToModify[rSublawIndex] == Constitutive::eInput::ENGINEERING_STRAIN)
            {
                Constitutive::eOutput derivativeOutputEnum =
                    GetDerivativeEnumSublaw(Constitutive::eOutput::ENGINEERING_STRAIN, itMainLawOutput.first);
                modifiedOutputMap.emplace(derivativeOutputEnum,
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(derivativeOutputEnum));
            }
            break;

        default:
            break;
        }
    }

    // Add sublaw outputs depending on main law inputs
    switch (mInputsToModify[rSublawIndex])
    {
    case Constitutive::eInput::ENGINEERING_STRAIN:
        for (const auto& itMainLawInput : rMainLawInputMap)
        {
            switch (itMainLawInput.first)
            {
            case Constitutive::eInput::ENGINEERING_STRAIN:
                modifiedOutputMap.emplace(Constitutive::eOutput::ENGINEERING_STRAIN,
                    ConstitutiveIOBase::makeConstitutiveIO<TDim>(Constitutive::eOutput::ENGINEERING_STRAIN));

                break;
            default:
                break;
            }
        }
        break;
    default:
        break;
    }
    return modifiedOutputMap;
}
