// $Id: ConstitutiveEngineeringStressStrain.cpp 112 2009-11-17 16:43:15Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/thermal/LinearHeatFlux.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient1D.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient2D.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux1D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux2D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::LinearHeatFlux::LinearHeatFlux() : ConstitutiveBase()
{
	mHeatCapacity = -1.;
    mThermalConductivity = -1.;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::LinearHeatFlux::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::LinearHeatFlux::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::LinearHeatFlux::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::LinearHeatFlux::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::LinearHeatFlux::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::LinearHeatFlux::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::LinearHeatFlux::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveHeatFlux" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveHeatFlux" << std::endl;
#endif
}

BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LinearHeatFlux)

#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::LinearHeatFlux::Evaluate3D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
	// check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    }

	// calculate engineering strain
	if(rConstitutiveInput.find(NuTo::Constitutive::eInput::TEMPERATURE_GRADIENT_3D)==rConstitutiveInput.end())
		throw MechanicsException("[NuTo::LinearHeatFlux::Evaluate] temperature gradient 3d needed.");
	const TemperatureGradient3D& temperatureGradient3D(rConstitutiveInput.find(NuTo::Constitutive::eInput::TEMPERATURE_GRADIENT_3D)->second->GetTemperatureGradient3D());

	//check, if an nonlinear iteration has to be performed, in this simple case, just calculate the linear elastic coefficients
    if (rConstitutiveOutput.find(NuTo::Constitutive::eOutput::HEAT_FLUX_3D)!=rConstitutiveOutput.end()
    		|| rConstitutiveOutput.find(NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D)!=rConstitutiveOutput.end())
    {
    	//in a nonlinear material routine, the Return mapping would be performed here
    }

    for (std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
    		itOutput != rConstitutiveOutput.end(); itOutput++)
    {
    	switch(itOutput->first)
    	{
    	case NuTo::Constitutive::eOutput::HEAT_FLUX_3D:
    	{
/*    		EngineeringStrain3D elasticEngineeringStrain(engineeringStrain);
    		// if temperature is an input, subtract thermal strains to get elastic strains
    		if (section->GetInputConstitutiveIsTemperature())
    		{
    			std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::eInput::TEMPERATURE));
    			if (itInput==rConstitutiveInput.end())
    				throw MechanicsException("[NuTo::LinearElasticEngineeringStress::Evaluate] temperature needed to evaluate thermal engineering strain3d.");
    			double temperature(itInput->second->GetTemperature());
    			double deltaStrain(mThermalExpansionCoefficient * temperature);
    			EngineeringStrain3D elasticEngineeringStrain;
    			elasticEngineeringStrain.mEngineeringStrain[0] -= deltaStrain;
    			elasticEngineeringStrain.mEngineeringStrain[1] -= deltaStrain;
    			elasticEngineeringStrain.mEngineeringStrain[2] -= deltaStrain;
    		}
			EngineeringStress3D& engineeringStress(itOutput->second->GetEngineeringStress3D());
		    // calculate Engineering stress
		    engineeringStress.mEngineeringStress[0] = C11 * elasticEngineeringStrain.mEngineeringStrain[0] + C12 * (elasticEngineeringStrain.mEngineeringStrain[1]+elasticEngineeringStrain.mEngineeringStrain[2]);
		    engineeringStress.mEngineeringStress[1] = C11 * elasticEngineeringStrain.mEngineeringStrain[1] + C12 * (elasticEngineeringStrain.mEngineeringStrain[0]+elasticEngineeringStrain.mEngineeringStrain[2]);
		    engineeringStress.mEngineeringStress[2] = C11 * elasticEngineeringStrain.mEngineeringStrain[2] + C12 * (elasticEngineeringStrain.mEngineeringStrain[0]+elasticEngineeringStrain.mEngineeringStrain[1]);
		    engineeringStress.mEngineeringStress[3] = C44 * elasticEngineeringStrain.mEngineeringStrain[3] ;
		    engineeringStress.mEngineeringStress[4] = C44 * elasticEngineeringStrain.mEngineeringStrain[4] ;
		    engineeringStress.mEngineeringStress[5] = C44 * elasticEngineeringStrain.mEngineeringStrain[5] ;
*/
    		break;

    	}
    	case NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D:
    	{
			ConstitutiveTangentLocal<3,3>& tangent(itOutput->second->AsConstitutiveTangentLocal_3x3());
 		    double *data(tangent.mTangent);

		    // store tangent at the output object
		    data[ 0] = mThermalConductivity;
		    data[ 1] = 0.;
		    data[ 2] = 0.;

		    data[ 6] = 0.;
		    data[ 7] = mThermalConductivity;
		    data[ 8] = 0.;

		    data[12] = 0.;
		    data[13] = 0.;
		    data[14] = mThermalConductivity;

		    tangent.SetSymmetry(true);
    		break;
    	}
    	default:
    		throw MechanicsException(std::string("[NuTo::LinearHeatFlux::Evaluate3D] output object)") +
    				NuTo::Constitutive::OutputToString(itOutput->first) +
    				std::string(" culd not be calculated, check the allocated material law and the section behavior."));
    	}
    }

    //update history variables but for linear , there is nothing to do

	return Error::SUCCESSFUL;
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::LinearHeatFlux::AllocateStaticDataEngineeringStress_EngineeringStrain3D(ElementBase* rElement)const
{
	return 0;
}

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::LinearHeatFlux::GetType() const
{
    return NuTo::Constitutive::LINEAR_HEAT_FLUX;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::LinearHeatFlux::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return true;
    case NuTo::Element::TETRAHEDRON4N:
        return true;
    case NuTo::Element::TETRAHEDRON10N:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if heat capacity is positive
//! @param rHeatCapacity ... heat capacity
void NuTo::LinearHeatFlux::CheckHeatCapacity(double rHeatCapacity) const
{
	if (rHeatCapacity < 0.0)
	{
		throw NuTo::MechanicsException("[NuTo::LinearHeatFlux::CheckDensity] The heat capacity must be a positive value.");
	}
}

//! @brief ... check if thermal conductivity is positive
//! @param rThermalConductivity ... thermal conductivity
void NuTo::LinearHeatFlux::CheckThermalConductivity(double rThermalConductivity) const
{
	if (rThermalConductivity < 0.0)
	{
		throw NuTo::MechanicsException("[NuTo::LinearHeatFlux::CheckDensity] The thermal conductivity must be a positive value.");
	}
}

// parameters /////////////////////////////////////////////////////////////

//! @brief ... get HeatCapacity
//! @return ... HeatCapacity
double NuTo::LinearHeatFlux::GetHeatCapacity() const
{
	return mHeatCapacity;
}


//! @brief ... set HeatCapacity
//! @param rE ... HeatCapacity
void NuTo::LinearHeatFlux::SetHeatCapacity(double rHeatCapacity)
{
    this->CheckHeatCapacity(rHeatCapacity);
    this->mHeatCapacity = rHeatCapacity;
    this->SetParametersValid();
}

//! @brief ... get ThermalConductivity
//! @return ... ThermalConductivity
double NuTo::LinearHeatFlux::GetThermalConductivity() const
{
	return this->mThermalConductivity;
}

//! @brief ... set density
//! @param rRho ... density
void NuTo::LinearHeatFlux::SetThermalConductivity(double rThermalConductivity)
{
    this->CheckThermalConductivity(rThermalConductivity);
    this->mThermalConductivity = rThermalConductivity;
    this->SetParametersValid();
}

// check parameters
void NuTo::LinearHeatFlux::CheckParameters()const
{
    this->CheckHeatCapacity(this->mHeatCapacity);
    this->CheckThermalConductivity(this->mThermalConductivity);
}

