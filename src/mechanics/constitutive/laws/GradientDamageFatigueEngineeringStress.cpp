#include <iostream>

#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "mechanics/constitutive/laws/GradientDamageFatigueEngineeringStress.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

double NuTo::GradientDamageFatigueEngineeringStress::GetParameterDouble(
        NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
    case Constitutive::eConstitutiveParameter::DENSITY:
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
    case Constitutive::eConstitutiveParameter::ENDURANCE_STRESS:
        return mEnduranceStress;
    case Constitutive::eConstitutiveParameter::FATIGUE_PARAMETER:
        return mFatigueParameter;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void NuTo::GradientDamageFatigueEngineeringStress::SetParameterDouble(
        NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
    case Constitutive::eConstitutiveParameter::DENSITY:
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
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
        NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData,
        const NuTo::ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput) const;
template std::pair<double, double> NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData<2>(
        NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData,
        const NuTo::ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput) const;
template std::pair<double, double> NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData<3>(
        NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData,
        const NuTo::ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput) const;


template <int TDim>
std::pair<double, double> NuTo::GradientDamageFatigueEngineeringStress::GetCurrentStaticData(
        NuTo::GradientDamageFatigueEngineeringStress::Data& rStaticData,
        const NuTo::ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput) const
{
    ////     1D completely explicit euler
    //    double kappa_old = rStaticData.GetData()[0];
    //
    //    double nonlocalEqStrain = rConstitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN)->operator[](0);
    //    double strain = rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->operator[](0);
    //    double sigmaEq = (1 - CalculateDamage(kappa_old))*mE*strain_old;
    //    double kappa = k(nonlocalEqStrain, sigmaEq, rStaticData.GetData());
    //
    //    if (rConstitutiveOutput.Contains(Constitutive::eOutput::UPDATE_STATIC_DATA))
    //        rStaticData.SetData(Eigen::Vector3d({kappa, nonlocalEqStrain, strain}));
    //
    //    double dKappa_dNonlocalEqStrain = dk_de(nonlocalEqStrain, sigmaEq, rStaticData.GetData());
    //
    //    return std::make_pair(kappa, dKappa_dNonlocalEqStrain);


    //    NuTo::ConstitutivePlaneState planeState(ePlaneState::PLANE_STRESS);
    //    if (TDim == 2)
    //        planeState =
    //        *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());

    double strain = rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->operator[](0);
    double nonlocalEqStrain = rConstitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN)->operator[](0);
    //    EngineeringStress<TDim> elasticStress =
    //    EngineeringStressHelper::GetStress(strain->AsEngineeringStrain<TDim>(), mE, mNu, planeState.GetPlaneState());

    double kappa = 0;
    double dKappa_dNonlocalEqStrain = 0;

    double kappa_old = rStaticData.GetData()[0];
    if (nonlocalEqStrain > kappa_old)
    {
        // static loading
        kappa = nonlocalEqStrain;
        dKappa_dNonlocalEqStrain = 1;
    }
    else
    {
        // other stuff
        kappa = -1; // damn.. just to enter the while
        double sigmaEq = 0;
        int iterations = 0;
        while (iterations < 1000 and std::abs(kappa_old - kappa) > 1.e-12)
        {
            ++iterations;
            kappa_old = kappa;
            sigmaEq = (1. - mDamageLaw->CalculateDamage(kappa)) * mE * strain;
            kappa = k(nonlocalEqStrain, sigmaEq, rStaticData.GetData());
        }
        if (iterations > 900)
        {
            std::cout << "Damn." << std::endl;
            throw;
        }

        double delta_e = nonlocalEqStrain - rStaticData.GetData()[1];
        if (delta_e < 0)
            dKappa_dNonlocalEqStrain = -F(sigmaEq);
        else
            dKappa_dNonlocalEqStrain = F(sigmaEq);
    }
    if (rConstitutiveOutput.Contains(Constitutive::eOutput::UPDATE_STATIC_DATA))
        rStaticData.SetData(Eigen::Vector2d({kappa, nonlocalEqStrain}));

    return std::make_pair(kappa, dKappa_dNonlocalEqStrain);
}

NuTo::Constitutive::eConstitutiveType NuTo::GradientDamageFatigueEngineeringStress::GetType() const
{
    return Constitutive::eConstitutiveType::GRADIENT_DAMAGE_FATIGUE_ENGINEERING_STRESS;
}
