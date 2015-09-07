#include "DryingShrinkage.h"

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include <nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h>
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFraction.h"

#include <math.h>



//! @brief Constructor
NuTo::DryingShrinkage::DryingShrinkage()
    : ConstitutiveBase()
{}



//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::DryingShrinkage::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::ELEMENT2D:
        return true;
    default:
        return false;
    }
}



//! @brief ... checks if the constitutive law has a specific parameter
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... true/false
bool NuTo::DryingShrinkage::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        default:
        {
            return false;
        }
    }
}



//! @brief ... checks if a constitutive law has an specific output
//! @return ... true/false
bool NuTo::DryingShrinkage::CheckOutputTypeCompatibility(NuTo::Constitutive::Output::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::Output::ENGINEERING_STRESS_2D_PORE_PRESSURE:
    case Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY_2D:
    case Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION_2D:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}



//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::DryingShrinkage::Evaluate2D(NuTo::ElementBase *rElement,
                                                      int rIp, const std::map<NuTo::Constitutive::Input::eInput,
                                                      const NuTo::ConstitutiveInputBase *> &rConstitutiveInput,
                                                      std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase *> &rConstitutiveOutput)
{
    // get section information determining which input on the constitutive level should be used
    //const SectionBase* section(rElement->GetSection());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }



    if(rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] water volume fraction needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] relative humidity needed to evaluate moisture transport.");
    }
    const RelativeHumidity&         relativeHumidity            (rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)              ->second->GetRelativeHumidity());
    const WaterVolumeFraction&      waterVolumeFraction         (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION)          ->second->GetWaterVolumeFraction());

    for (auto itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch(itOutput->first)
        {

        case NuTo::Constitutive::Output::ENGINEERING_STRESS_2D_PORE_PRESSURE:
        {
            EngineeringStress2D& engineeringStressPorePressure(itOutput->second->GetEngineeringStress2D());

            engineeringStressPorePressure[0]    = (//mAtmosphericPressure
                                                  - waterVolumeFraction(0) / mPorosity *
                                                  mDensityWater * mIdealGasConstant * mTemperature / mMolarMassWater *
                                                  log(relativeHumidity(0)));
            engineeringStressPorePressure[1]    = engineeringStressPorePressure[0];
            engineeringStressPorePressure[2]    = 0.0;
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY_2D:
        {
            ConstitutiveTangentLocal<2,1>&  tangent_D_EngineeringStress_D_RH(itOutput->second->AsConstitutiveTangentLocal_2x1());
            tangent_D_EngineeringStress_D_RH[0] = - waterVolumeFraction(0) / mPorosity *
                                                  mDensityWater * mIdealGasConstant * mTemperature /
                                                  (mMolarMassWater * relativeHumidity(0));
            tangent_D_EngineeringStress_D_RH[1] = tangent_D_EngineeringStress_D_RH[0];
            break;
        }
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION_2D:
        {
            ConstitutiveTangentLocal<2,1>&  tangent_D_EngineeringStress_D_WV(itOutput->second->AsConstitutiveTangentLocal_2x1());
            tangent_D_EngineeringStress_D_WV[0] = - mDensityWater * mIdealGasConstant * mTemperature /
                                                  (mMolarMassWater * mPorosity) *
                                                  log(relativeHumidity(0));
            tangent_D_EngineeringStress_D_WV[1] = tangent_D_EngineeringStress_D_WV[0];
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



//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::DryingShrinkage::GetType() const
{
    return Constitutive::eConstitutiveType::DRYING_SHRINKAGE;
}



//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::DryingShrinkage::HaveTmpStaticData() const
{
    return false;
}



//! @brief ... check parameters of the constitutive relationship
//! if one check fails, an exception is thrown
void NuTo::DryingShrinkage::CheckParameters() const
{

}
