#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/laws/HeatConduction.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"

#include "nuto/mechanics/elements/ElementBase.h"

NuTo::HeatConduction::HeatConduction() : ConstitutiveBase()
{
    mK = 0.0;
    mCt = 0.0;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
template<class Archive>
void NuTo::HeatConduction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize HeatConduction" << std::endl;
#endif
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
//    & BOOST_SERIALIZATION_NVP(mE)
//    & BOOST_SERIALIZATION_NVP(mNu)
//    & BOOST_SERIALIZATION_NVP(mRho)
//    & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize HeatConduction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::HeatConduction)
#endif // ENABLE_SERIALIZATION

NuTo::ConstitutiveInputMap NuTo::HeatConduction::GetConstitutiveInputs(
    const ConstitutiveOutputMap& rConstitutiveOutput,
    const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::HEAT_FLUX:
        {
            constitutiveInputMap[Constitutive::Input::TEMPERATURE_GRADIENT];
            break;
        }
        case NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        {
            // for nonlinear:
            //constitutiveInputMap[Constitutive::Input::TEMPERATURE];
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, Constitutive::OutputToString(itOutput.first) + " cannot be calculated by this constitutive law.");
        }
    }

    return constitutiveInputMap;
}

NuTo::Error::eError NuTo::HeatConduction::Evaluate1D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::HEAT_FLUX:
        {
            const auto& temperatureGradient = *rConstitutiveInput.at(Constitutive::Input::TEMPERATURE_GRADIENT);
            ConstitutiveIOBase& heatFlux = *itOutput.second;
            heatFlux[0] = mK * temperatureGradient[0];
            break;
        }
        case NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1,1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0,0) = mK;
            break;
        }
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            //nothing to be done for update routine
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Output object "
                    + Constitutive::OutputToString(itOutput.first)
                    + " could not be calculated, check the allocated material law and the section behavior.");
        }
    }
    return Error::SUCCESSFUL;
}

NuTo::Error::eError NuTo::HeatConduction::Evaluate2D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    return Error::SUCCESSFUL;
}

NuTo::Error::eError NuTo::HeatConduction::Evaluate3D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    return Error::SUCCESSFUL;
}

bool NuTo::HeatConduction::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
        case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
        {
            return true;
        }
        default:
        {
            return false;
        }
    }
}

double NuTo::HeatConduction::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
        return this->mK;
    case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
        return this->mCt;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void NuTo::HeatConduction::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    ConstitutiveBase::CheckParameterDouble(rIdentifier, rValue);
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
        this->mK = rValue;
        break;
    case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
        this->mCt = rValue;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

bool NuTo::HeatConduction::CheckOutputTypeCompatibility(NuTo::Constitutive::Output::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::Output::HEAT_FLUX:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}

NuTo::Constitutive::eConstitutiveType NuTo::HeatConduction::GetType() const
{
    return NuTo::Constitutive::HEAT_CONDUCTION;
}

bool NuTo::HeatConduction::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::CONTINUUMELEMENT:
        return true;
    default:
        return false;
    }
}

void NuTo::HeatConduction::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Thermal conductivity          : " << this->mK << "\n";
    rLogger << "    Heat capacity                 : " << this->mCt << "\n";
}

void NuTo::HeatConduction::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, mK);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::HEAT_CAPACITY, mCt);
}
