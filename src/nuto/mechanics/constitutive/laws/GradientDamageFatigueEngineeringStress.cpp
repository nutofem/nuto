#include "nuto/mechanics/constitutive/laws/GradientDamageFatigueEngineeringStress.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

double NuTo::GradientDamageFatigueEngineeringStress::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        case Constitutive::eConstitutiveParameter::DENSITY:
        case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
        case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
        case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
            return GradientDamageEngineeringStress::GetParameterDouble(rIdentifier);
        case Constitutive::eConstitutiveParameter::ENDURANCE_STRESS:
            return mEnduranceStress;
        case Constitutive::eConstitutiveParameter::FATIGUE_PARAMETER:
            return mFatigueParameter;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__,"Constitutive law does not have the requested variable");
    }
}

void NuTo::GradientDamageFatigueEngineeringStress::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        case Constitutive::eConstitutiveParameter::DENSITY:
        case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
        case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
        case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
            GradientDamageEngineeringStress::SetParameterDouble(rIdentifier, rValue);
        case Constitutive::eConstitutiveParameter::ENDURANCE_STRESS:
            mEnduranceStress = rValue;
            break;
        case Constitutive::eConstitutiveParameter::FATIGUE_PARAMETER:
            mFatigueParameter = rValue;
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
    SetParametersValid();
}

double NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData(NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData,
                                                                          const NuTo::ConstitutiveInputMap& rConstitutiveInput) const
{
    return 0;
}

NuTo::Constitutive::eConstitutiveType NuTo::GradientDamageFatigueEngineeringStress::GetType() const
{
    return Constitutive::eConstitutiveType::GRADIENT_DAMAGE_FATIGUE_ENGINEERING_STRESS;
}
