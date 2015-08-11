// $Id$
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
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

// constructor
NuTo::ConstitutiveBase::ConstitutiveBase()
{
    this->mParametersValid = false;
}

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::ConstitutiveBase::Evaluate1D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput)
{
	std::cout << "[ConstitutiveBase::Evaluate1D] make this function pure virtual." << std::endl;
	throw std::string("[ConstitutiveBase::Evaluate] make this function pure virtual.");
	return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::ConstitutiveBase::Evaluate2D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput)
{
	std::cout << "[ConstitutiveBase::Evaluate1D] make this function pure virtual." << std::endl;
	throw std::string("[ConstitutiveBase::Evaluate] make this function pure virtual.");
	return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::ConstitutiveBase::Evaluate3D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput)
{
	std::cout << "[ConstitutiveBase::Evaluate3D] make this function pure virtual." << std::endl;
	throw std::string("[ConstitutiveBase::Evaluate] make this function pure virtual.");
    return Error::SUCCESSFUL;
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
bool NuTo::ConstitutiveBase::GetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterBool] This constitutive law has no variables of type bool.");
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterBool] This constitutive law has no variables of type bool.");
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
double NuTo::ConstitutiveBase::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterDouble] This constitutive law has no variables of type double.");
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetParameterDouble] This constitutive law has no variables of type double.");
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::ConstitutiveBase::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetParameterFullVectorDouble] This constitutive law has no variables of type NuTo::FullVector<double, Eigen::Dynamic>.");
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveBase::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::SetParameterFullVectorDouble] This constitutive law has no variables of type NuTo::FullVector<double, Eigen::Dynamic>.");
}


//! @brief ... get factor to modify Young's modulus (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorYoungsModulus(const ElementBase* rElement,int rIp) const
{
    return 1;
}


//! @brief ... get factor to modify Poisson's ratio (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorPoissonsRatio(const ElementBase* rElement,int rIp) const
{
    return 1;
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

//! @brief ... get factor to modify yield strength (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorYieldStrength(const ElementBase* rElement,int rIp) const
{
    return 1;
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

//! @brief ... get factor to modify hardening modulus (using random fields)
//! @param rElement ...  element
//! @param rIp ...  integration point
double NuTo::ConstitutiveBase::GetRanfieldFactorHardeningModulus(const ElementBase* rElement,int rIp) const
{
    return 1;
}



//! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
//! @param rRelativeHumidity ... relative humidity
//! @return ... equilibrium water volume fraction
double NuTo::ConstitutiveBase::GetEquilibriumWaterVolumeFraction(double rRelativeHumidity, NuTo::FullVector<double,Eigen::Dynamic> rCoeffs) const
{
    throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::GetEquilibriumWaterVolumeFraction] The constitutive relationship does not have this parameter.");
}




// modify parameter validity flag
void NuTo::ConstitutiveBase::SetParametersValid()
{
    try
    {
        this->CheckParameters();
    }
    catch (NuTo::MechanicsException)
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

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveBase::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement)const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Allocate routine for 1D EngineeringStressStrain not implemented.");
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveBase::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement)const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AllocateStaticDataEngineeringStress_EngineeringStrain2D] Allocate routine for 2D EngineeringStressStrain not implemented.");
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::ConstitutiveBase::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement)const
{
	throw NuTo::MechanicsException("[NuTo::ConstitutiveBase::AllocateStaticDataEngineeringStress_EngineeringStrain3D] Allocate routine for 3D EngineeringStressStrain not implemented.");
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
