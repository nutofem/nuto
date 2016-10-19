// $Id$
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/math/FullMatrix.h"
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

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief ... checks if the constitutive law has a specific parameter
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... true/false
bool NuTo::ConstitutiveBase::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::CheckHaveParameter] Not implemented for this constitutive law.");
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
bool NuTo::ConstitutiveBase::GetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterBool] This constitutive law has no variables of type bool.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterBool] This constitutive law has no variables of type bool.");
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::ConstitutiveBase::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterDouble] This constitutive law has no variables of type double.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetParameterDouble] This constitutive law has no variables of type double.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterFunction(std::function<std::array<double, 2>(double)>)
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "This constitutive law has no variables of type double.");
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::ConstitutiveBase::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterFullVectorDouble] This constitutive law has no variables of type NuTo::FullVector<double, Eigen::Dynamic>.");
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetParameterFullVectorDouble] This constitutive law has no variables of type NuTo::FullVector<double, Eigen::Dynamic>.");
}

//! @brief checks parameters, throws if the check failed
void NuTo::ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
    case Constitutive::eConstitutiveParameter::DENSITY:
    case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
    case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
    {
        if (rValue < 0.)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, Constitutive::ConstitutiveParameterToString(rIdentifier) + " must be > 0. (value: " + std::to_string(rValue) + ").");
    }
    break;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
    {
        if (rValue <= -1.0)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Poisson's ratio must be greater or equal to -1.0 (value: " + std::to_string(rValue) + ").");
        if (rValue >= 0.5)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Poisson's ratio must be smaller or equal to 0.5 (value: " + std::to_string(rValue) + ").");
    }
    break;

    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
    break;

    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "material parameter check not implemented.");
    }
}

//! @brief ... get yield strength for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding yield strength
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::ConstitutiveBase::GetYieldStrength() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... add yield strength
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  yield strength
void NuTo::ConstitutiveBase::AddYieldStrength(double rEpsilon, double rSigma)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AddYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... get hardening modulus for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding hardening modulus
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::ConstitutiveBase::GetHardeningModulus() const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}

//! @brief ... add hardening modulus
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  hardening modulus
void NuTo::ConstitutiveBase::AddHardeningModulus(double rEpsilon, double rH)
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AddYieldStrength] The constitutive relationship does not have a parameter yield strength.");
}


//! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
//! @param rRelativeHumidity ... relative humidity
//! @return ... equilibrium water volume fraction
double NuTo::ConstitutiveBase::GetEquilibriumWaterVolumeFraction(double rRelativeHumidity, NuTo::FullVector<double,Eigen::Dynamic> rCoeffs) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetEquilibriumWaterVolumeFraction] The constitutive relationship does not have this parameter.");
}



//! @brief ... checks if a constitutive law has an specific output
//! @return ... true/false
bool NuTo::ConstitutiveBase::CheckOutputTypeCompatibility(Constitutive::eOutput rOutputEnum) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::CheckOutputTypeCompatibility] Function not implemented for this constitutive law.");
}

// modify parameter validity flag
void NuTo::ConstitutiveBase::SetParametersValid()
{
    try
    {
        this->CheckParameters();
    }
    catch (NuTo::MechanicsException& e)
    {
        this->mParametersValid = false;
        return;
    }
    catch (...)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetParametersValid] Unhandled exception");
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
