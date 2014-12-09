#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterPhaseFraction.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#include "nuto/mechanics/constitutive/moistureTransport/MoistureTransport.h"

NuTo::MoistureTransport::MoistureTransport()
    :   ConstitutiveBase()
{
    SetParametersValid();
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool                                        NuTo::MoistureTransport::CheckElementCompatibility          (Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return false;
    case NuTo::Element::PLANE2D10N:
        return false;
    case NuTo::Element::PLANE2D15N:
        return false;
    case NuTo::Element::PLANE2D3N:
        return false;
    case NuTo::Element::PLANE2D4N:
        return false;
    case NuTo::Element::PLANE2D4NSPECTRALORDER2:
        return false;
    case NuTo::Element::PLANE2D4NSPECTRALORDER3:
        return false;
    case NuTo::Element::PLANE2D4NSPECTRALORDER4:
        return false;
    case NuTo::Element::PLANE2D6N:
        return false;
    case NuTo::Element::TETRAHEDRON4N:
        return false;
    case NuTo::Element::TETRAHEDRON10N:
        return false;
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return false;
    default:
        return false;
    }
}

//! @brief ... check if the mass exchange rate is non-negative
//! @param rMassExchangeRate ... mass exchange rate
void                                        NuTo::MoistureTransport::CheckMassExchangeRate              (double rMassExchangeRate) const
{
    if (rMassExchangeRate < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckMassExchangeRate] The mass exchange rate must have a non-negative value.");
    }
}

//! @brief ... check parameters of the constitutive relationship
//! if one check fails, an exception is thrwon
void                                        NuTo::MoistureTransport::CheckParameters                    () const
{
    CheckMassExchangeRate(mR);
    CheckPorosity(mEpsP);
    CheckVaporPhaseDiffusionCoefficient(mDV);
    CheckVaporPhaseDiffusionExponent(mAlphaV);
    CheckVaporPhaseSaturationDensity(mRhoVS);
    CheckWaterPhaseDensity(mRhoW);
    CheckWaterPhaseDiffusionCoefficient(mDW);
    CheckWaterPhaseDiffusionExponent(mAlphaW);
}

//! @brief ... check if the porosity is a value between 0 and 1
//! @param rPorosity ... porosity
void                                        NuTo::MoistureTransport::CheckPorosity                      (double rPorosity) const
{
    if (rPorosity <= 0 || rPorosity >= 1)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckPorosity] The porosity must have a value between 0 and 1.");
    }
}

//! @brief ... check if vapor phase diffusion coefficient is non-negative
//! @param rVaporPhaseDiffusionCoefficient ... vapor phase diffusion coefficient
void                                        NuTo::MoistureTransport::CheckVaporPhaseDiffusionCoefficient(double rVaporPhaseDiffusionCoefficient) const
{
    if (rVaporPhaseDiffusionCoefficient < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckVaporPhaseDiffusionCoefficient] The vapor phase diffusion coefficient must have a non-negative value.");
    }
}

//! @brief ... check if vapor phase diffusion exponent is non-negative
//! @param rVaporPhaseDiffusionExponent ... vapor phase diffusion exponent
void                                        NuTo::MoistureTransport::CheckVaporPhaseDiffusionExponent   (double rVaporPhaseDiffusionExponent) const
{
    if (rVaporPhaseDiffusionExponent < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckVaporPhaseDiffusionExponent] The vapor phase diffusion exponent must have a non-negative value.");
    }
}

//! @brief ... check if vapor phase saturation density is positive and non-zero
//! @param rVaporPhaseSaturationDensity ... vapor phase saturation density
void                                        NuTo::MoistureTransport::CheckVaporPhaseSaturationDensity   (double rVaporPhaseSaturationDensity) const
{
    if (rVaporPhaseSaturationDensity <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckVaporPhaseSaturationDensity] The vapor phase saturation density must be a non-negative, non-zero value.");
    }
}

//! @brief ... check if water phase density is positive and non-zero
//! @param rWaterPhaseDensity ... water phase density
void                                        NuTo::MoistureTransport::CheckWaterPhaseDensity             (double rWaterPhaseDensity) const
{
    if (rWaterPhaseDensity <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckWaterPhaseDensity] The water phase density must be a non-negative, non-zero value.");
    }
}

//! @brief ... check if water phase diffusion coefficient is non-negative
//! @param rWaterPhaseDiffusionCoefficient ... water phase diffusion coefficient
void                                        NuTo::MoistureTransport::CheckWaterPhaseDiffusionCoefficient(double rWaterPhaseDiffusionCoefficient) const
{
    if (rWaterPhaseDiffusionCoefficient < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckWaterPhaseDiffusionCoefficient] The water phase diffusion coefficient must have a non-negative value.");
    }
}

//! @brief ... check if water phase diffusion exponent is non-negative
//! @param rWaterPhaseDiffusionExponent ... water phase diffusion exponent
void                                        NuTo::MoistureTransport::CheckWaterPhaseDiffusionExponent   (double rWaterPhaseDiffusionExponent) const
{
    if (rWaterPhaseDiffusionExponent < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckWaterPhaseDiffusionExponent] The water phase diffusion exponent must have a non-negative value.");
    }
}

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError                         NuTo::MoistureTransport::Evaluate1D                         (ElementBase* rElement, int rIp,
                                                                                                         const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                                         std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    // get section information determining which input on the constitutive level should be used
    const SectionBase* section(rElement->GetSection());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }



    if(rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_PHASE_FRACTION)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] water phase fraction needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] relative humidity needed to evaluate moisture transport.");
    }

    const RelativeHumidity&     relativeHumidity    (rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)->second->GetRelativeHumidity());
    const WaterPhaseFraction&   waterPhaseFraction  (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_PHASE_FRACTION)->second->GetWaterPhaseFraction());

    for (auto itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch(itOutput->first)
        {
        case NuTo::Constitutive::Output::VAPOR_PHASE_DIFFUSION_COEFFICIENT:
        {
            ConstitutiveTangentLocal<1,1>& VaporPhaseDiffusionCoefficient(itOutput->second->AsConstitutiveTangentLocal_1x1());
            VaporPhaseDiffusionCoefficient(0,0)=mDV * pow(1 - (waterPhaseFraction(0) / mEpsP), mAlphaV);
            VaporPhaseDiffusionCoefficient.SetSymmetry(true);
            break;
        }
        case NuTo::Constitutive::Output::WATER_PHASE_DIFFUSION_COEFFICIENT:
        {
            ConstitutiveTangentLocal<1,1>& WaterPhaseDiffusionCoefficient(itOutput->second->AsConstitutiveTangentLocal_1x1());
            WaterPhaseDiffusionCoefficient(0,0)=mDW * pow(waterPhaseFraction(0) / mEpsP, mAlphaW);
            WaterPhaseDiffusionCoefficient.SetSymmetry(true);
            break;
        }
        case NuTo::Constitutive::Output::PHASE_MASS_EXCHANGE_RATE:
        {
            ConstitutiveTangentLocal<1,1>& PhaseMassExchangeRate(itOutput->second->AsConstitutiveTangentLocal_1x1());
            PhaseMassExchangeRate(0,0)=mR;
            PhaseMassExchangeRate.SetSymmetry(true);
            break;
        }
        case NuTo::Constitutive::Output::PHASE_MASS_EXCHANGE_RATE_TIMES_EQUILIBRIUM_SORPTION_CURVE:
        {
            ConstitutiveTangentLocal<1,1>& PhaseMassExchangeRateTimesEquilibriumSorptionCurve(itOutput->second->AsConstitutiveTangentLocal_1x1());
            PhaseMassExchangeRateTimesEquilibriumSorptionCurve(0,0) = mR * (mActualSorptionCoeff(0) +
                                                                            mActualSorptionCoeff(1) * relativeHumidity(0) +
                                                                            mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0));
            PhaseMassExchangeRateTimesEquilibriumSorptionCurve.SetSymmetry(true);
            break;
        }
        case NuTo::Constitutive::Output::WATER_PHASE_DENSITY:
        {
            ConstitutiveTangentLocal<1,1>& WaterPhaseDensity(itOutput->second->AsConstitutiveTangentLocal_1x1());
            WaterPhaseDensity(0,0) = mRhoW;
            WaterPhaseDensity.SetSymmetry(true);
            break;
        }
        case NuTo::Constitutive::Output::VAPOR_PHASE_SATURATION_DENSITY_TIMES_VAPOR_PHASE_VOLUME_FRACTION:
        {
            if (waterPhaseFraction(0)>mEpsP)
            {
                throw MechanicsException(std::string("[NuTo::MoistureTransport::Evaluate1D] Water phase volume fraction bigger than porosity!"));
            }
            ConstitutiveTangentLocal<1,1>& VaporPhaseSaturationDensityTimesVaporPhaseVolumeFraction(itOutput->second->AsConstitutiveTangentLocal_1x1());
            VaporPhaseSaturationDensityTimesVaporPhaseVolumeFraction(0,0) = mRhoVS * (mEpsP  - waterPhaseFraction(0));
            VaporPhaseSaturationDensityTimesVaporPhaseVolumeFraction.SetSymmetry(true);
            break;
        }
        case NuTo::Constitutive::Output::VAPOR_PHASE_SATURATION_DENSITY_TIMES_RELATIVE_HUMIDITY:
        {
            ConstitutiveTangentLocal<1,1>& VaporPhaseSaturationDensityTimesRelativeHumidity(itOutput->second->AsConstitutiveTangentLocal_1x1());
            VaporPhaseSaturationDensityTimesRelativeHumidity(0,0) = mRhoVS * relativeHumidity(0);
            VaporPhaseSaturationDensityTimesRelativeHumidity.SetSymmetry(true);
            break;
        }
        default:
            throw MechanicsException(std::string("[NuTo::MoistureTransport::Evaluate1D ] output object)") +
                                     NuTo::Constitutive::OutputToString(itOutput->first) +
                                     std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    return Error::SUCCESSFUL;
}

//! @brief ... get mass exchange rate between water phase and vapor phase
//! @return ... mass exchange rate
double                                      NuTo::MoistureTransport::GetMassExchangeRate                () const
{
    return mR;
}

//! @brief ... get porosity
//! @return ... porosity
double                                      NuTo::MoistureTransport::GetPorosity                        () const
{
    return mEpsP;
}

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType       NuTo::MoistureTransport::GetType                            () const
{
    return NuTo::Constitutive::MOISTURE_TRANSPORT;
}

//! @brief ... get vapor phase diffusion coefficient
//! @return ... vapor phase diffusion coefficient
double                                      NuTo::MoistureTransport::GetVaporPhaseDiffusionCoefficient  () const
{
    return mDV;
}

//! @brief ... get vapor phase diffusion exponent
//! @return ... vapor phase diffusion exponent
double                                      NuTo::MoistureTransport::GetVaporPhaseDiffusionExponent     () const
{
    return mAlphaV;
}

//! @brief ... get vapor phase saturation density
//! @return ... vapor phase saturation density
double                                      NuTo::MoistureTransport::GetVaporPhaseSaturationDensity     () const
{
    return mRhoVS;
}

//! @brief ... get water phase density
//! @return ... water phase density
double                                      NuTo::MoistureTransport::GetWaterPhaseDensity               () const
{
    return mRhoW;
}

//! @brief ... get water phase diffusion coefficient
//! @return ... water phase diffusion coefficient
double                                      NuTo::MoistureTransport::GetWaterPhaseDiffusionCoefficient  () const
{
    return mDW;
}

//! @brief ... get water phase diffusion exponent
//! @return ... water phase diffusion exponent
double                                      NuTo::MoistureTransport::GetWaterPhaseDiffusionExponent     () const
{
    return mAlphaW;
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool                                        NuTo::MoistureTransport::HaveTmpStaticData                  () const
{
    return false;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void                                        NuTo::MoistureTransport::Info                               (unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Mass exchange rate                 : " << this->GetMassExchangeRate()               << "\n";
    rLogger << "    Porosity                           : " << this->GetPorosity()                       << "\n";
    rLogger << "    Vapor phase diffusion coefficient  : " << this->GetVaporPhaseDiffusionCoefficient() << "\n";
    rLogger << "    Vapor phase diffusion exponent     : " << this->GetVaporPhaseDiffusionExponent()    << "\n";
    rLogger << "    Vapor phase saturation density     : " << this->GetVaporPhaseSaturationDensity()    << "\n";
    rLogger << "    Water phase density                : " << this->GetWaterPhaseDensity()              << "\n";
    rLogger << "    Water phase diffusion coefficient  : " << this->GetWaterPhaseDiffusionCoefficient() << "\n";
    rLogger << "    Water phase diffusion exponent     : " << this->GetWaterPhaseDiffusionExponent()    << "\n";
}

//! @brief ... set mass exchange rate between water phase and vapor phase
//! @param ... mass exchange rate
void                                        NuTo::MoistureTransport::SetMassExchangeRate                (double rMassExchangeRate)
{
    CheckMassExchangeRate(rMassExchangeRate);
    mR = rMassExchangeRate;
    SetParametersValid();
}

//! @brief ... set porosity
//! @param ... porosity
void                                        NuTo::MoistureTransport::SetPorosity                        (double rPorosity)
{
    CheckPorosity(rPorosity);
    mEpsP = rPorosity;
    SetParametersValid();
}

//! @brief ... set vapor phase diffusion coefficient
//! @param ... vapor phase diffusion coefficient
void                                        NuTo::MoistureTransport::SetVaporPhaseDiffusionCoefficient  (double rVaporPhaseDiffusionCoefficient)
{
    CheckVaporPhaseDiffusionCoefficient(rVaporPhaseDiffusionCoefficient);
    mDV = rVaporPhaseDiffusionCoefficient;
    SetParametersValid();
}

//! @brief ... set vapor phase diffusion exponent
//! @param ... vapor phase diffusion exponent
void                                        NuTo::MoistureTransport::SetVaporPhaseDiffusionExponent     (double rVaporPhaseDiffusionExponent)
{
    CheckVaporPhaseDiffusionExponent(rVaporPhaseDiffusionExponent);
    mAlphaV = rVaporPhaseDiffusionExponent;
    SetParametersValid();
}

//! @brief ... set vapor phase saturation density
//! @param ... vapor phase saturation density
void                                        NuTo::MoistureTransport::SetVaporPhaseSaturationDensity     (double rVaporPhaseSaturationDensity)
{
    CheckVaporPhaseSaturationDensity(rVaporPhaseSaturationDensity);
    mRhoVS = rVaporPhaseSaturationDensity;
    SetParametersValid();
}

//! @brief ... set water phase density
//! @param ... water phase density
void                                        NuTo::MoistureTransport::SetWaterPhaseDensity               (double rWaterPhaseDensity)
{
    CheckWaterPhaseDensity(rWaterPhaseDensity);
    mRhoW = rWaterPhaseDensity;
    SetParametersValid();
}

//! @brief ... set water phase diffusion coefficient
//! @param ... water phase diffusion coefficient
void                                        NuTo::MoistureTransport::SetWaterPhaseDiffusionCoefficient  (double rWaterPhaseDiffusionCoefficient)
{
    CheckWaterPhaseDiffusionCoefficient(rWaterPhaseDiffusionCoefficient);
    mDW = rWaterPhaseDiffusionCoefficient;
    SetParametersValid();
}

//! @brief ... set water phase diffusion exponent
//! @param ... water phase diffusion exponent
void                                        NuTo::MoistureTransport::SetWaterPhaseDiffusionExponent     (double rWaterPhaseDiffusionExponent)
{
    CheckWaterPhaseDiffusionExponent(rWaterPhaseDiffusionExponent);
    mAlphaW = rWaterPhaseDiffusionExponent;
    SetParametersValid();
}



