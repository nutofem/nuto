#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

using namespace NuTo::Constitutive;

void NuTo::AdditiveInputExplicit::AddConstitutiveLaw(NuTo::ConstitutiveBase& rConstitutiveLaw,
                                                     Constitutive::eInput rModiesInput)
{
    if (rModiesInput == Constitutive::eInput::NONE)
    {
        if (mMainLaw != nullptr)
            throw MechanicsException(
                    __PRETTY_FUNCTION__,
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


NuTo::ConstitutiveInputMap
NuTo::AdditiveInputExplicit::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap& rConstitutiveOutput) const
{
    // Get Inputs for output returning constitutive law
    ConstitutiveInputMap mainLawConstitutiveInputMap(
            mMainLaw->GetConstitutiveInputs(rConstitutiveOutput));

    // Get Inputs for input modifying constitutive laws
    ConstitutiveInputMap sublawsConstitutiveInputMap;
    for (size_t i = 0; i < mSublaws.size(); ++i)
    {

        // Get the necessary inputs for the sublaws ---> the modifications to the main laws inputs are returned as
        // output from the sublaws. Therefore one first needs the sublaws output list, depending on the main laws
        // inputs and the global outputs.

        // INFO regarding the template function call: The dimension can be chosen freely, because the created object is
        // never used. Should be 1 because of the lowest construction costs.  Alternative would be to copy the called
        // function and replace all objects with nullptr. But then you have to maintain 2 nearly identical functions.
        ConstitutiveOutputMap sublawOutputMap =
                GetSublawOutputMap<1>(mainLawConstitutiveInputMap, rConstitutiveOutput, i);

        // Don't merge the sublaw inputs directly into the main laws input map!  ---> When more than one sublaw is
        // attached the inputs of the first might effect the following laws outputs - have a look at the line above!
        sublawsConstitutiveInputMap.Merge(mSublaws[i]->GetConstitutiveInputs(sublawOutputMap));
    }
    return mainLawConstitutiveInputMap.Merge(sublawsConstitutiveInputMap);
}


NuTo::Constitutive::eConstitutiveType NuTo::AdditiveInputExplicit::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT;
}


NuTo::Constitutive::eOutput
NuTo::AdditiveInputExplicit::GetDerivativeEnumSublaw(NuTo::Constitutive::eOutput rParameter,
                                                     NuTo::Constitutive::eOutput rMainDerivative) const
{
    switch (rParameter)
    {
    case Constitutive::eOutput::ENGINEERING_STRAIN:
        switch (rMainDerivative)
        {
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
            return Constitutive::eOutput::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY;

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
            return Constitutive::eOutput::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION;

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
            return Constitutive::eOutput::D_STRAIN_D_TEMPERATURE;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "No partial derivative defined for parameter " +
                                                                  Constitutive::OutputToString(rParameter) +
                                                                  " and global derivative " +
                                                                  Constitutive::OutputToString(rMainDerivative));
        }
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "No partial derivatives defined for parameter " +
                                                              Constitutive::OutputToString(rParameter));
    }
}


template <int TDim>
NuTo::ConstitutiveOutputMap
NuTo::AdditiveInputExplicit::GetSublawOutputMap(const NuTo::ConstitutiveInputMap& rMainLawInputMap,
                                                const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                                int rSublawIndex) const
{
    ConstitutiveOutputMap modifiedOutputMap;

    // Add sublaw outputs depending on main law outputs
    for (const auto& itMainLawOutput : rMainLawOutputMap)
    {
        switch (itMainLawOutput.first)
        {
        case Constitutive::eOutput::SHRINKAGE_STRAIN_VISUALIZE:
        case Constitutive::eOutput::THERMAL_STRAIN:
            modifiedOutputMap.emplace(itMainLawOutput.first, itMainLawOutput.second->clone());
            break;
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
            if (mInputsToModify[rSublawIndex] == Constitutive::eInput::ENGINEERING_STRAIN)
            {
                Constitutive::eOutput derivativeOutputEnum =
                        GetDerivativeEnumSublaw(Constitutive::eOutput::ENGINEERING_STRAIN, itMainLawOutput.first);
                modifiedOutputMap.emplace(derivativeOutputEnum,
                                          ConstitutiveIOBase::makeConstitutiveIO<TDim>(derivativeOutputEnum));
            }
            break;
        case Constitutive::eOutput::UPDATE_STATIC_DATA:
            modifiedOutputMap[Constitutive::eOutput::UPDATE_STATIC_DATA];
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
                modifiedOutputMap.emplace(
                        Constitutive::eOutput::ENGINEERING_STRAIN,
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

template NuTo::ConstitutiveOutputMap
NuTo::AdditiveInputExplicit::GetSublawOutputMap<1>(const NuTo::ConstitutiveInputMap& rMainLawInputMap,
                                                   const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                                   int rSublawIndex) const;
template NuTo::ConstitutiveOutputMap
NuTo::AdditiveInputExplicit::GetSublawOutputMap<2>(const NuTo::ConstitutiveInputMap& rMainLawInputMap,
                                                   const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                                   int rSublawIndex) const;
template NuTo::ConstitutiveOutputMap
NuTo::AdditiveInputExplicit::GetSublawOutputMap<3>(const NuTo::ConstitutiveInputMap& rMainLawInputMap,
                                                   const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                                   int rSublawIndex) const;
