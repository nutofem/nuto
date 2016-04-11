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
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
    & BOOST_SERIALIZATION_NVP(mK)
    & BOOST_SERIALIZATION_NVP(mCt);
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
            throw MechanicsException(__PRETTY_FUNCTION__, Constitutive::OutputToString(itOutput.first)
                    + " cannot be calculated by this constitutive law.");
        }
    }

    return constitutiveInputMap;
}

template<int TDim>
NuTo::Error::eError NuTo::HeatConduction::Evaluate(NuTo::ElementBase *rElement,
        int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    auto eye = Eigen::MatrixXd::Identity(TDim, TDim);

    InputData<TDim> inputData;
    for (auto itInput : rConstitutiveInput)
    {
        switch(itInput.first)
        {
        case NuTo::Constitutive::Input::TEMPERATURE_GRADIENT:
            inputData.mTemperatureGradient = (*static_cast<ConstitutiveVector<TDim>*>(itInput.second)).AsVector();
            break;
        case NuTo::Constitutive::Input::CALCULATE_STATIC_DATA:
        case NuTo::Constitutive::Input::TIME_STEP:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Input object "
                + Constitutive::InputToString(itInput.first) + " is not needed by the heat conduction law.");
        }
    }

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::HEAT_FLUX:
        {
            Eigen::Matrix<double, TDim, 1>& heatFlux = (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second));
            heatFlux = mK * eye * inputData.mTemperatureGradient;
            break;
        }
        case NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        {
            Eigen::Matrix<double, TDim, TDim>& tangent =
                (*static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second));
            tangent = mK * eye;
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
