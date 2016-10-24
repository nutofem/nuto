#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/EngineeringStressHelper.h"
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


template std::pair<double, double> NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData<1>(
    NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData, const NuTo::ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput) const;
template std::pair<double, double> NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData<2>(
    NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData, const NuTo::ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput) const;
template std::pair<double, double> NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData<3>(
    NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData, const NuTo::ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput) const;


template <int TDim>
std::pair<double, double> NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData(
    NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData,
    const NuTo::ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput) const
{
    NuTo::ConstitutivePlaneState planeState(ePlaneState::PLANE_STRESS);
    if (TDim == 2)
        planeState = *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());

    auto& strain = rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN);
    EngineeringStress<TDim> engineeringStress = EngineeringStressHelper::GetStress(strain->AsEngineeringStrain<TDim>(), mE, mNu, planeState.GetPlaneState());
    double kappa_old = rStaticData.GetData()[0];
    engineeringStress *= (1 - CalculateDamage(kappa_old));
    double sigmaEq = engineeringStress.GetVonMisesStress(planeState.GetPlaneState());

    // This is explicit euler.
    double nonlocalEqStrain = rConstitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN)->operator[](0);
    double kappa = k(nonlocalEqStrain, sigmaEq, rStaticData.GetData());
    double dKappa_dNonlocalEqStrain = dk_de(nonlocalEqStrain, sigmaEq, rStaticData.GetData());

    if (rConstitutiveOutput.Contains(Constitutive::eOutput::UPDATE_STATIC_DATA))
        rStaticData.SetData(Eigen::Vector2d({kappa, nonlocalEqStrain}));

    return std::make_pair(kappa, dKappa_dNonlocalEqStrain);
}

NuTo::Constitutive::eConstitutiveType NuTo::GradientDamageFatigueEngineeringStress::GetType() const
{
    return Constitutive::eConstitutiveType::GRADIENT_DAMAGE_FATIGUE_ENGINEERING_STRESS;
}
