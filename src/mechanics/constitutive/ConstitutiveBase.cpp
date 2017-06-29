// $Id$
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "base/Logger.h"
#include "base/Exception.h"

//! @brief ... checks if the constitutive law has a specific parameter
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... true/false
bool NuTo::ConstitutiveBase::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::CheckHaveParameter] Not implemented for this constitutive law.");
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
bool NuTo::ConstitutiveBase::GetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::GetParameterBool] This constitutive law has no variables of type bool.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::GetParameterBool] This constitutive law has no variables of type bool.");
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::ConstitutiveBase::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::GetParameterDouble] This constitutive law has no variables of type double.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::SetParameterDouble] This constitutive law has no variables of type double.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterFunction(std::function<std::array<double, 2>(double)>)
{
    throw NuTo::Exception(__PRETTY_FUNCTION__, "This constitutive law has no variables of type double.");
}

void NuTo::ConstitutiveBase::SetDamageLaw(std::shared_ptr<NuTo::Constitutive::DamageLaw> damageLaw)
{
    throw NuTo::Exception(__PRETTY_FUNCTION__, "This constitutive law has no damage law.");
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
Eigen::VectorXd NuTo::ConstitutiveBase::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::GetParameterFullVectorDouble] This constitutive law has no variables of type Eigen::VectorXd.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, Eigen::VectorXd rValue)
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::SetParameterFullVectorDouble] This constitutive law has no variables of type Eigen::VectorXd.");
}

//! @brief checks parameters, throws if the check failed
void NuTo::ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
    case Constitutive::eConstitutiveParameter::DENSITY:
    case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
    case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
    {
        if (rValue < 0.)
            throw NuTo::Exception(__PRETTY_FUNCTION__, Constitutive::ConstitutiveParameterToString(rIdentifier) + " must be > 0. (value: " + std::to_string(rValue) + ").");
    }
    break;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
    {
        if (rValue <= -1.0)
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Poisson's ratio must be greater or equal to -1.0 (value: " + std::to_string(rValue) + ").");
        if (rValue >= 0.5)
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Poisson's ratio must be smaller or equal to 0.5 (value: " + std::to_string(rValue) + ").");
    }
    break;

    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
    break;

    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "material parameter check not implemented.");
    }
}


//! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
//! @param rRelativeHumidity ... relative humidity
//! @return ... equilibrium water volume fraction
double NuTo::ConstitutiveBase::GetEquilibriumWaterVolumeFraction(double rRelativeHumidity, Eigen::VectorXd rCoeffs) const
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::GetEquilibriumWaterVolumeFraction] The constitutive relationship does not have this parameter.");
}



//! @brief ... checks if a constitutive law has an specific output
//! @return ... true/false
bool NuTo::ConstitutiveBase::CheckOutputTypeCompatibility(Constitutive::eOutput rOutputEnum) const
{
    throw NuTo::Exception("[NuTo::ConstitutiveBase::CheckOutputTypeCompatibility] Function not implemented for this constitutive law.");
}

// modify parameter validity flag
void NuTo::ConstitutiveBase::SetParametersValid()
{
    try
    {
        this->CheckParameters();
    }
    catch (NuTo::Exception& e)
    {
        this->mParametersValid = false;
        return;
    }
    catch (...)
    {
        throw NuTo::Exception("[NuTo::ConstitutiveBase::SetParametersValid] Unhandled exception");
    }
    this->mParametersValid = true;
}


// info routine
void NuTo::ConstitutiveBase::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    std::cout << "    parameter validity flag: " << this->mParametersValid << std::endl;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize constitutive Base" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mParametersValid);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Constitutive Base" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstitutiveBase)
#endif // ENABLE_SERIALIZATION
