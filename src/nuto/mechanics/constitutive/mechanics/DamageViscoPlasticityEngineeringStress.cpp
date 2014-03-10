// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/Logger.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataDamageViscoPlasticity3D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/DamageViscoPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

#include <math.h>

#define sqrt3 1.732050808

NuTo::DamageViscoPlasticityEngineeringStress::DamageViscoPlasticityEngineeringStress() : ConstitutiveBase()
{
	mE = 0.;
	mNu = 0.;
	mRho = 0.;
	mTensileStrength = 0.;
	mCompressiveStrength = 0.;
	mBiaxialCompressiveStrength = 0.;
	mViscosity = 1.;
	mDamageDistribution = 1.;
	mViscoplasticYieldSurfaceOffset = 0.;
    mFractureEnergy = 0.;
    mDamage = true;
	SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::DamageViscoPlasticityEngineeringStress::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
   std::cout << "start serialize DamageViscoPlasticityEngineeringStress" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
      & BOOST_SERIALIZATION_NVP(mE)
      & BOOST_SERIALIZATION_NVP(mNu)
      & BOOST_SERIALIZATION_NVP(mRho)
      & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient)
      & BOOST_SERIALIZATION_NVP(mTensileStrength)
      & BOOST_SERIALIZATION_NVP(mCompressiveStrength)
      & BOOST_SERIALIZATION_NVP(mBiaxialCompressiveStrength)
      & BOOST_SERIALIZATION_NVP(mViscosity)
      & BOOST_SERIALIZATION_NVP(mDamageDistribution)
      & BOOST_SERIALIZATION_NVP(mViscoplasticYieldSurfaceOffset)
      & BOOST_SERIALIZATION_NVP(mFractureEnergy);
#ifdef DEBUG_SERIALIZATION
   std::cout << "finish serialize DamageViscoPlasticityEngineeringStress" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::DamageViscoPlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::DamageViscoPlasticityEngineeringStress::Evaluate1D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate1D] not implemented for 1D.");
}
// Case of Elasticity 1D;
//{
//	// get section information determining which input on the constitutive level should be used
//	const SectionBase* section(rElement->GetSection());
//
//	// check if parameters are valid
//    if (this->mParametersValid == false)
//    {
//   		//throw an exception giving information related to the wrong parameter
//    	CheckParameters();
//    }
//
//	EngineeringStrain1D engineeringStrain;
//	// calculate engineering strain
//	if(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D)==rConstitutiveInput.end())
//		throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate] deformation gradient 1d needed to evaluate engineering strain2d.");
//	const DeformationGradient1D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D)->second->GetDeformationGradient1D());
//	deformationGradient.GetEngineeringStrain(engineeringStrain);
//
//    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
//    		itOutput != rConstitutiveOutput.end(); itOutput++)
//    {
//    	switch(itOutput->first)
//    	{
//    	case NuTo::Constitutive::Output::ENGINEERING_STRESS_1D:
//    	{
//    		EngineeringStrain1D elasticEngineeringStrain(engineeringStrain);
//    		// if temperature is an input, subtract thermal strains to get elastic strains
//    		if (section->GetInputConstitutiveIsTemperature())
//    		{
//    			std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
//    			if (itInput==rConstitutiveInput.end())
//    				throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate1D] temperature needed to evaluate thermal engineering strain1d.");
//    			double temperature(itInput->second->GetTemperature());
//    			double deltaStrain(mThermalExpansionCoefficient * temperature);
//    			EngineeringStrain1D elasticEngineeringStrain;
//    			elasticEngineeringStrain[0] -= deltaStrain;
//    		}
//			EngineeringStress1D& engineeringStress(itOutput->second->GetEngineeringStress1D());
//			// calculate Engineering stress
//			engineeringStress = mE*elasticEngineeringStrain;
//
//		    break;
//    	}
//    	case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
//    	{
//    		//this is for the visualize routines
//    		EngineeringStrain1D elasticEngineeringStrain(engineeringStrain);
//    		// if temperature is an input, subtract thermal strains to get elastic strains
//    		if (section->GetInputConstitutiveIsTemperature())
//    		{
//    			std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
//    			if (itInput==rConstitutiveInput.end())
//    				throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate2D] temperature needed to evaluate thermal engineering strain2d.");
//    			double temperature(itInput->second->GetTemperature());
//    			double deltaStrain(mThermalExpansionCoefficient * temperature);
//    			EngineeringStrain1D elasticEngineeringStrain;
//    			elasticEngineeringStrain[0] -= deltaStrain;
//    		}
//			EngineeringStress3D& engineeringStress(itOutput->second->GetEngineeringStress3D());
//
//			// calculate Engineering stress
//			engineeringStress[0] = mE * elasticEngineeringStrain[0];
//			engineeringStress[1] = 0.;
//			engineeringStress[2] = 0.;
//			engineeringStress[3] = 0.;
//			engineeringStress[4] = 0.;
//			engineeringStress[5] = 0.;
//
//		    break;
//    	}
//    	case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D:
//    	{
//			ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
// 		    tangent(0,0)=mE;
//		    tangent.SetSymmetry(true);
//    		break;
//    	}
//    	case NuTo::Constitutive::Output::ENGINEERING_STRAIN_1D:
//    	{
//    		EngineeringStrain1D& engineeringStrain1D(itOutput->second->GetEngineeringStrain1D());
//			engineeringStrain1D[0] = engineeringStrain[0];
//    	}
//    	break;
//    	case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
//    	{
//    		EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
//			engineeringStrain3D[0] = engineeringStrain[0];
//			engineeringStrain3D[1] = -mNu*engineeringStrain[0];
//			engineeringStrain3D[2] = engineeringStrain3D[1];
//			engineeringStrain3D[3] = 0.;
//			engineeringStrain3D[4] = 0.;
//			engineeringStrain3D[5] = 0.;
//    	}
//    	break;
//    	case NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D:
//    	{
//    		EngineeringStrain3D& engineeringPlasticStrain(itOutput->second->GetEngineeringStrain3D());
//    		engineeringPlasticStrain[0] = 0.;
//    		engineeringPlasticStrain[1] = 0.;
//    		engineeringPlasticStrain[2] = 0.;
//    		engineeringPlasticStrain[3] = 0.;
//    		engineeringPlasticStrain[4] = 0.;
//    		engineeringPlasticStrain[5] = 0.;
//    		break;
//    	}
//    	case NuTo::Constitutive::Output::DAMAGE:
//    	{
//    		itOutput->second->GetDamage().SetDamage(0.);
//    		break;
//    	}
//    	case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
//    	case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
//    	{
//    	    //nothing to be done for update routine
//    		break;
//    	}
//    	default:
//    		throw MechanicsException(std::string("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate1D] output object)") +
//    				NuTo::Constitutive::OutputToString(itOutput->first) +
//    				std::string(" could not be calculated, check the allocated material law and the section behavior."));
//    	}
//    }
//
//    //update history variables but for linear elastic, there is nothing to do
//
//	return Error::SUCCESSFUL;
//}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::DamageViscoPlasticityEngineeringStress::Evaluate2D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate2D] not implemented for 2D.");
}
// Case of Elasticity 2D;
//{
//	// get section information determining which input on the constitutive level should be used
//	const SectionBase* section(rElement->GetSection());
//
//	// check if parameters are valid
//    if (this->mParametersValid == false)
//    {
//   		//throw an exception giving information related to the wrong parameter
//    	CheckParameters();
//    }
//
//	EngineeringStrain2D engineeringStrain;
//	// calculate engineering strain
//	if(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D)==rConstitutiveInput.end())
//		throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate] deformation gradient 2d needed to evaluate engineering strain2d.");
//	const DeformationGradient2D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D)->second->GetDeformationGradient2D());
//	deformationGradient.GetEngineeringStrain(engineeringStrain);
//
//    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
//    		itOutput != rConstitutiveOutput.end(); itOutput++)
//    {
//    	switch(itOutput->first)
//    	{
//    	case NuTo::Constitutive::Output::ENGINEERING_STRESS_2D:
//    	{
//    		EngineeringStrain2D elasticEngineeringStrain(engineeringStrain);
//    		// if temperature is an input, subtract thermal strains to get elastic strains
//    		if (section->GetInputConstitutiveIsTemperature())
//    		{
//    			std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
//    			if (itInput==rConstitutiveInput.end())
//    				throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate2D] temperature needed to evaluate thermal engineering strain2d.");
//    			double temperature(itInput->second->GetTemperature());
//    			double deltaStrain(mThermalExpansionCoefficient * temperature);
//    			EngineeringStrain2D elasticEngineeringStrain;
//    			elasticEngineeringStrain[0] -= deltaStrain;
//    			elasticEngineeringStrain[1] -= deltaStrain;
//    		}
//			EngineeringStress2D& engineeringStress(itOutput->second->GetEngineeringStress2D());
//		    // calculate Engineering stress
//
//		    switch(rElement->GetSection()->GetType())
//		    {
//		    case Section::PLANE_STRAIN:{
//				// calculate coefficients of the material matrix
//				double C11, C12, C33;
//				this->CalculateCoefficients3D(C11, C12, C33);
//
//				// calculate Engineering stress
//				engineeringStress[0] = C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1];
//				engineeringStress[1] = C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0];
//				engineeringStress[2] = C33 * elasticEngineeringStrain[2] ;
//		    	break;}
//		    case Section::PLANE_STRESS:{
//				// calculate coefficients of the material matrix
//				double C11, C12, C33;
//				this->CalculateCoefficients2DPlainStress(C11, C12, C33);
//
//				// calculate Engineering stress
//				engineeringStress[0] = C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1];
//				engineeringStress[1] = C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0];
//				engineeringStress[2] = C33 * elasticEngineeringStrain[2] ;
//		    	break;}
//		    default:
//		    	throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Invalid type of 2D section behavoir found!!!");
//		    }
//
//
//		    break;
//    	}
//    	case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
//    	{
//    		//this is for the visualize routines
//    		EngineeringStrain2D elasticEngineeringStrain(engineeringStrain);
//    		// if temperature is an input, subtract thermal strains to get elastic strains
//    		if (section->GetInputConstitutiveIsTemperature())
//    		{
//    			std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
//    			if (itInput==rConstitutiveInput.end())
//    				throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate2D] temperature needed to evaluate thermal engineering strain2d.");
//    			double temperature(itInput->second->GetTemperature());
//    			double deltaStrain(mThermalExpansionCoefficient * temperature);
//    			EngineeringStrain2D elasticEngineeringStrain;
//    			elasticEngineeringStrain[0] -= deltaStrain;
//    			elasticEngineeringStrain[1] -= deltaStrain;
//    		}
//			EngineeringStress3D& engineeringStress(itOutput->second->GetEngineeringStress3D());
//		    // calculate Engineering stress
//
//		    switch(rElement->GetSection()->GetType())
//		    {
//		    case Section::PLANE_STRAIN:{
//				// calculate coefficients of the material matrix
//				double C11, C12, C33;
//				this->CalculateCoefficients3D(C11, C12, C33);
//
//				// calculate Engineering stress
//				engineeringStress[0] = C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1];
//				engineeringStress[1] = C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0];
//				engineeringStress[2] = C12 * (elasticEngineeringStrain[0]+elasticEngineeringStrain[1]);
//				engineeringStress[3] = C33 * elasticEngineeringStrain[2] ;
//				engineeringStress[4] = 0.;
//				engineeringStress[5] = 0.;
//		    	break;}
//		    case Section::PLANE_STRESS:{
//				// calculate coefficients of the material matrix
//				double C11, C12, C33;
//				this->CalculateCoefficients2DPlainStress(C11, C12, C33);
//
//				// calculate Engineering stress
//				engineeringStress[0] = C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1];
//				engineeringStress[1] = C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0];
//				engineeringStress[2] = 0.;
//				engineeringStress[3] = C33 * elasticEngineeringStrain[2] ;
//				engineeringStress[4] = 0.;
//				engineeringStress[5] = 0.;
//		    	break;}
//		    default:
//		    	throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStressFromEngineeringStrain] Invalid type of 2D section behavoir found!!!");
//		    }
//
//
//		    break;
//    	}
//    	case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D:
//    	{
//			ConstitutiveTangentLocal<3,3>& tangent(itOutput->second->AsConstitutiveTangentLocal_3x3());
//
// 			assert(rElement->GetSection()!=0);
// 		    switch(rElement->GetSection()->GetType())
// 		    {
// 		    case Section::PLANE_STRAIN:{
// 			    // calculate coefficients of the material matrix
// 			    double C11, C12, C33;
// 			    this->CalculateCoefficients3D(C11, C12, C33);
//
// 			    // store tangent at the output object
// 			    tangent(0,0) = C11;
// 			    tangent(1,0) = C12;
// 			    tangent(2,0) = 0;
//
// 			    tangent(0,1) = C12;
// 			    tangent(1,1) = C11;
// 			    tangent(2,1) = 0;
//
// 			    tangent(0,2) = 0.;
// 			    tangent(1,2) = 0.;
// 			    tangent(2,2) = C33;
// 		    	break;}
// 		    case Section::PLANE_STRESS:{
// 			    // calculate coefficients of the material matrix
// 			    double C11, C12, C33;
// 				this->CalculateCoefficients2DPlainStress(C11, C12, C33);
//
// 			    // store tangent at the output object
// 			    tangent(0,0) = C11;
// 			    tangent(1,0) = C12;
// 			    tangent(2,0) = 0;
//
// 			    tangent(0,1) = C12;
// 			    tangent(1,1) = C11;
// 			    tangent(2,1) = 0;
//
// 			    tangent(0,2) = 0.;
// 			    tangent(1,2) = 0.;
// 			    tangent(2,2) = C33;
// 		   	break;}
// 		    default:
// 		    	throw MechanicsException("[NuTo::LinearElastic::GetTangent_EngineeringStress_EngineeringStrain] Invalid type of 2D section behavoir found!!!");
// 		    }
//
//		    tangent.SetSymmetry(true);
//    		break;
//    	}
//    	case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
//    	{
//    		EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
//    	    switch(rElement->GetSection()->GetType())
//    	    {
//    	    case Section::PLANE_STRAIN:
//    	    	engineeringStrain3D[0] = engineeringStrain[0];
//    	    	engineeringStrain3D[1] = engineeringStrain[1];
//    	    	engineeringStrain3D[2] = 0;
//    	    	engineeringStrain3D[3] = engineeringStrain[2];
//    	    	engineeringStrain3D[4] = 0.;
//    	    	engineeringStrain3D[5] = 0.;
//    	    	break;
//    	    case Section::PLANE_STRESS:
//    	    	engineeringStrain3D[0] = engineeringStrain[0];
//    	    	engineeringStrain3D[1] = engineeringStrain[1];
//    	    	engineeringStrain3D[2] = mNu/(mNu-1.)*(engineeringStrain[0]+engineeringStrain[1]);
//    	    	engineeringStrain3D[3] = engineeringStrain[2];
//    	    	engineeringStrain3D[4] = 0.;
//    	    	engineeringStrain3D[5] = 0.;
//    	    	break;
//    	    default:
//    	    	throw MechanicsException("[NuTo::LinearElastic::GetEngineeringStrainFromEngineeringStrain] Invalid type of 2D section behavoir found!!!");
//    	    }
//    	}
//    	break;
//    	case NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D:
//    	{
//    		EngineeringStrain3D& engineeringPlasticStrain(itOutput->second->GetEngineeringStrain3D());
//    		engineeringPlasticStrain[0] = 0.;
//    		engineeringPlasticStrain[1] = 0.;
//    		engineeringPlasticStrain[2] = 0.;
//    		engineeringPlasticStrain[3] = 0.;
//    		engineeringPlasticStrain[4] = 0.;
//    		engineeringPlasticStrain[5] = 0.;
//    		break;
//    	}
//    	case NuTo::Constitutive::Output::DAMAGE:
//    	{
//    		itOutput->second->GetDamage().SetDamage(0.);
//    		break;
//    	}
//    	case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
//    	case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
//    	{
//    	    //nothing to be done for update routine
//    		break;
//    	}
//    	default:
//    		throw MechanicsException(std::string("[NuTo::DamageViscoPlasticityEngineeringStressEngineeringStress::Evaluate3D] output object)") +
//    				NuTo::Constitutive::OutputToString(itOutput->first) +
//    				std::string(" culd not be calculated, check the allocated material law and the section behavior."));
//    	}
//    }
//
//    //update history variables but for linear elastic, there is nothing to do
//
//	return Error::SUCCESSFUL;
//}

//! @brief ... evaluate the constitutive relation in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::DamageViscoPlasticityEngineeringStress::Evaluate3D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
// THIS IS TEST NEWTON
//	NuTo::FullVector<double,Eigen::Dynamic> x(3);
//	bool check;
//	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> (DamageViscoPlasticityEngineeringStress::*fdjacAn)
//			(NuTo::FullVector<double,Eigen::Dynamic>) const;
//	fdjacAn = &DamageViscoPlasticityEngineeringStress::DResidualAn;
//	x[0] = -1., x[1] = 2., x[2] = 3.;
//	std::cout << x.transpose() << std::endl;
//	this->Newton(x, check, &DamageViscoPlasticityEngineeringStress::DResidualAn);
//	this->Newton(x, check);
//	std::cout << x.transpose() << std::endl;
//	std::cout << check << std::endl;
// THIS IS TEST NEWTON

	// get section information determining which input on the constitutive level should be used
	const SectionBase* section(rElement->GetSection());

	// check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    }

	EngineeringStrain3D engineeringStrain;

	// calculate engineering strain via Deformation gradient
//	if(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D)==rConstitutiveInput.end())
//		throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate] deformation gradient 3d needed to evaluate engineering strain3d.");
//	const DeformationGradient3D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D)->second->GetDeformationGradient3D());
//	deformationGradient.GetEngineeringStrain(engineeringStrain);

	// get engineering strain directly
	if(rConstitutiveInput.find(NuTo::Constitutive::Input::ENGINEERING_STRAIN_3D)==rConstitutiveInput.end())
		throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate] engineering strain 3d needed.");
	engineeringStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::ENGINEERING_STRAIN_3D)->second->GetEngineeringStrain3D();

	//subtract thermal strain
	EngineeringStrain3D MechanicEngineeringStrain(engineeringStrain);
//    cout << "Eval.: Def Grad" << endl;
//    deformationGradient.Info();
    cout << "Eval.: eing Strain" << engineeringStrain.transpose() << endl;
	// if temperature is an input, subtract thermal strains to get elastic strains
	if (section->GetInputConstitutiveIsTemperature())
	{
		std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
		if (itInput==rConstitutiveInput.end())
			throw MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate3D] temperature needed to evaluate thermal engineering strain3d.");
		double temperature(itInput->second->GetTemperature());
		double deltaStrain(mThermalExpansionCoefficient * temperature);
		MechanicEngineeringStrain[0] -= deltaStrain;
		MechanicEngineeringStrain[1] -= deltaStrain;
		MechanicEngineeringStrain[2] -= deltaStrain;
	}
    cout << "Eval.: Mech Eng Strain" << MechanicEngineeringStrain.transpose() << endl;

	// declare output of Return Mapping: stress, algorithmic tangent and static data
	EngineeringStress3D engineeringStress;
	ConstitutiveTangentLocal<6,6>* tangent(0);
	ConstitutiveStaticDataDamageViscoPlasticity3D newStaticData;

	// check whether Return Mapping should be done: check output
	bool performReturnMapping(false);

	auto itStress = rConstitutiveOutput.find(NuTo::Constitutive::Output::ENGINEERING_STRESS_3D);
	if (itStress!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}
	auto itStiffness = rConstitutiveOutput.find(NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D);
	if (itStiffness!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
		tangent = &(itStiffness->second->AsConstitutiveTangentLocal_6x6());
	}
	auto itUpdate = rConstitutiveOutput.find(NuTo::Constitutive::Output::UPDATE_STATIC_DATA);
	if (itUpdate!=rConstitutiveOutput.end())
	{
		performReturnMapping = true;
	}

	// check, if an nonlinear iteration has to be performed
    if (performReturnMapping)
    {
    	// perform return mapping
    	NuTo::Error::eError errorReturnMapping = ReturnMapping3D(rElement, rIp, MechanicEngineeringStrain, &engineeringStress, tangent, &newStaticData ,rElement->GetStructure()->GetLogger());
    	if (errorReturnMapping!=Error::SUCCESSFUL)
    		return errorReturnMapping;
    }

    // calculate output
    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
    		itOutput != rConstitutiveOutput.end(); itOutput++)
    {
    	switch(itOutput->first)
    	{
    	case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
    	{
    		if (mDamage == true) {
    			itOutput->second->GetEngineeringStress3D() = (1. - newStaticData.mOmegaCompr)*engineeringStress;
			} else {
				itOutput->second->GetEngineeringStress3D() = engineeringStress;
			}
   		break;
    	}
    	case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D:
    	{
    		//tangent has already been calculated in Return Mapping
			itOutput->second->AsConstitutiveTangentLocal_6x6().SetSymmetry(true);

    		break;
    	}
    	case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
    		itOutput->second->GetEngineeringStrain3D() = engineeringStrain;
    		break;
    	case NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D:
    	{
    		EngineeringStrain3D& engineeringPlasticStrain(itOutput->second->GetEngineeringStrain3D());
    		engineeringPlasticStrain = newStaticData.mEpsilonP;
    		break;
    	}
    	case NuTo::Constitutive::Output::ENGINEERING_VISCOPLASTIC_STRAIN_3D:
    	{
    		EngineeringStrain3D& engineeringViscoPlasticStrain(itOutput->second->GetEngineeringStrain3D());
    		engineeringViscoPlasticStrain = newStaticData.mEpsilonVp;
    		break;
    	}
    	case NuTo::Constitutive::Output::DAMAGE:
    	{
    		itOutput->second->GetDamage().SetDamage(newStaticData.mOmegaCompr);
   		break;
    	}
    	case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
    	break;
    	case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
    	{
    	    *(rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3D()) = newStaticData;
    	    rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3D()->mPrevStrain = MechanicEngineeringStrain;
    	    rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3D()->mPrevSigma = engineeringStress;
    	}
		break;
    	default:
    		throw MechanicsException(std::string("[NuTo::DamageViscoPlasticityEngineeringStress::Evaluate3D] output object)") +
    				NuTo::Constitutive::OutputToString(itOutput->first) +
    				std::string(" culd not be calculated, check the allocated material law and the section behavior."));
    	}
    }


	return Error::SUCCESSFUL;
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::DamageViscoPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement)const
{
	return 0;
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::DamageViscoPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement)const
{
	return 0;
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::DamageViscoPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement)const
{
	return new ConstitutiveStaticDataDamageViscoPlasticity3D();
}

//! @brief ... performs the return mapping procedure in 3D
//! @param rElement ... structure
//! @param rIp ... integration point
//! @param rEngineeringStrain ... mechanical strain at the end of time increment
//! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
//! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
//! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
#define sqrt_2div3 0.81649658
#define toleranceYieldSurface mTOLF/100.  //tolerance whether a point is on the yield surface or not

NuTo::Error::eError NuTo::DamageViscoPlasticityEngineeringStress::ReturnMapping3D(const ElementBase* rElement,int rIp,
		const EngineeringStrain3D& rEngineeringStrain,
		EngineeringStress3D* rNewStress,
		ConstitutiveTangentLocal<6,6>* rNewTangent,
		ConstitutiveStaticDataDamageViscoPlasticity3D* rNewStaticData,
		Logger& rLogger) const
{
	// get structure information determining the beginning and the end of the time increment
	const StructureBase* structure(rElement->GetStructure());

	// times
	double InitTime, EndTime, DeltaTime;
	InitTime = structure->GetPrevTime();
	EndTime  = structure->GetTime();
	DeltaTime = EndTime - InitTime;
	cout << "*** InitTime = " << InitTime << "EndTime = " << EndTime << "Time Increment = " << DeltaTime << endl;

	// get material parameters
    double f_ct  = mTensileStrength;
    double f_c1  = mCompressiveStrength;
    double f_c2  = mBiaxialCompressiveStrength;
    double Hoffset = mViscoplasticYieldSurfaceOffset;
    double viscosity = mViscosity;
    double alpha = mDamageDistribution;

    assert(f_c2 > f_c1);
    assert(f_c1 > 0.);
    assert(f_c2 > 0.);
    assert(Hoffset <= 0.);
    assert(viscosity > 0.);
    assert(alpha >= 0. && alpha <= 1.);

    double Beta = sqrt3*(f_c2-f_c1) / (2*f_c2-f_c1);
    double HP  = f_c2*f_c1 / (sqrt3*(2*f_c2-f_c1));     // H for the plastic yield surface
    double HVP = HP + Hoffset;							// H for the viscoplastic yield surface

    // get elastic matrix
	// calculate coefficients of the linear elastic material matrix
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ElasticStiffness(6, 6);
	double C11, C12, C44;
	this->CalculateCoefficients3D(C11, C12, C44);
	ElasticStiffness << C11, C12, C12,  0., 0.,  0.,
						C12, C11, C12,  0., 0.,  0.,
						C12, C12, C11,  0., 0.,  0.,
						 0.,  0.,  0., C44, 0.,  0.,
						 0.,  0.,  0.,  0., C44, 0.,
						 0.,  0.,  0.,  0., 0., C44;

    // get state variables to begin of the time increment
    const ConstitutiveStaticDataDamageViscoPlasticity3D *OldStaticData = rElement->GetStaticData(rIp)->AsDamageViscoPlasticity3D();
    // extract total mechanical strain to begin of the time increment
    const EngineeringStress3D prevStress(OldStaticData->GetPrevStress());
    const EngineeringStrain3D prevStrain(OldStaticData->GetPrevStrain());  // prevStrain is the mechanical strain, without thermal strain

    // calculate trial stress at the end of time increment
    EngineeringStress3D TrialStress;
    TrialStress = ElasticStiffness*(rEngineeringStrain - OldStaticData->mEpsilonP - OldStaticData->mEpsilonVp);
    cout << "rEing Strain" << rEngineeringStrain.transpose() << endl;
    cout << "Old plastic strain = " << (OldStaticData->mEpsilonP).transpose() <<
    		"Old viscoplastic strain = " <<(OldStaticData->mEpsilonVp).transpose()<< endl;
    cout << "trial stress" << TrialStress.transpose() << endl;

    // calculate elastic solution
    if (rNewStress!=0) {
    	(*rNewStress) = TrialStress;
    }
    if (rNewTangent!=0) {
    	(*rNewTangent) = ElasticStiffness;
    }
    if (rNewStaticData!=0) {
    	(*rNewStaticData) = (*OldStaticData);
    	rNewStaticData->mVP = TrialStress.YieldSurfaceDruckerPrager3D(Beta, HP)/(DeltaTime*this->mE);
    }
    // if inelastic case, that is the viscoplastic yield surface criterion provides a positive value
    if (TrialStress.YieldSurfaceDruckerPrager3D(Beta, HVP) > -toleranceYieldSurface) {

    	// compound the vector of unknowns
    	NuTo::FullVector<double,Eigen::Dynamic> Unknown(14);

    	// initialize start values
    	Unknown.segment<6>(0) = ElasticStiffness*(rEngineeringStrain - prevStrain); 	// stress increment = Unknown(0:5)
    	Unknown.segment<6>(6) = ElasticStiffness*(rEngineeringStrain - prevStrain); 	// stress projection increment = Unknown(6:11)

    	Unknown[12] = TrialStress.YieldSurfaceDruckerPrager3D(Beta, HP)/(DeltaTime*this->mE);	// mVP = Unknown(12) plastic state variable
    	Unknown[13] = TrialStress.YieldSurfaceDruckerPrager3D(Beta, HVP)/(DeltaTime*this->mE);	// mVviscoP = Unknown(13) viscoplastic state variable

    	// compose vector of known Parameter, which are necessary for ResidualAn
    	NuTo::FullVector<double,Eigen::Dynamic> Parameter(17);

    	for (int i = 0; i < 6; i++) {
    		Parameter[i] = rEngineeringStrain[i] - prevStrain[i]; 	// Parameter(0:5) increment of mechanical strain
    		Parameter[i+6] = prevStress[i];   					   	// Parameter(6:11) stress at the beginning of the time increment
		}

    	Parameter[12] = DeltaTime;	                 		       // Parameter(12) time increment itself
    	Parameter[13] = Beta;									   // Parameter(13) Beta
    	Parameter[14] = HP;										   // Paramater(14) HP
    	Parameter[15] = HVP;									   // Paramater(15) HVP
    	Parameter[16] = viscosity;								   // Paramater(16) viscosity

    	const NuTo::FullVector<double,Eigen::Dynamic> ParameterList(Parameter);

        // prepare starting Newton solver with respect to the "Unknown"
        bool check;
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> (DamageViscoPlasticityEngineeringStress::*fdjacAn)
        	(const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>) const;

        // set Jacobi to analytical Jacobi
        fdjacAn = &DamageViscoPlasticityEngineeringStress::DResidualAn;

        // start Newton solver
        try
        {
        this->Newton(ParameterList,Unknown,check,fdjacAn);
//    	this->Newton(ParameterList,Unknown,check);   // use for numerical Jacobi
//        cout << "Unknown from Newton" << Unknown.transpose() << endl;;
//        cout << "Parameter from Newton" << ParameterList.transpose()<< endl;
		std::cout << boolalpha << std::endl;
        cout << "check from Newton" << check;

        }
        catch (...)
        {
            rLogger << "[NuTo::DamageViscoPlasticityEngineeringStress::ReturnMapping3D] No convergence after MAXITS steps, check the Newton." << "\n";

            return Error::NO_CONVERGENCE;
        }

        // calculate inelastic solution
		EngineeringStress3D NewStressTemp, NewProjectionStress;
    	for (int i = 0; i < 6; i++) {
    		NewStressTemp[i] = prevStress[i] + Unknown[i];
    		NewProjectionStress[i] = prevStress[i] + Unknown[i+6];
		}
		NuTo::FullMatrix<double,6,6> rd2F_d2Sigma;
		NuTo::FullVector<double,6> rdF_dSigma;

        NewStressTemp.YieldSurfaceDruckerPrager3DDerivatives(rdF_dSigma,rd2F_d2Sigma, Beta);

        // calculate stress
        if (rNewStress!=0) {
        	(*rNewStress) = NewStressTemp;
        	(*rNewStress).Info(1);
        }

        // calculate static data
        if (rNewStaticData!=0) {
        	EngineeringStrain3D DeltaEpsilonP, DeltaEpsilonVp;			// increments of plastic and viscoplastic strain respectively
        	//plastic state variable
        	rNewStaticData->mVP = Unknown[12];							// plastic state variable
        	// plastic strain
        	DeltaEpsilonP = Plus(rNewStaticData->mVP)*DeltaTime*rdF_dSigma;
        	rNewStaticData->mEpsilonP = OldStaticData->mEpsilonP + DeltaEpsilonP;
            //viscoplastic strain
            double YieldViscoPlastic = NewStressTemp.YieldSurfaceDruckerPrager3D(Beta, HVP);
            DeltaEpsilonVp = (ElasticStiffness.fullPivLu().solve(NewStressTemp - NewProjectionStress)).eval();
            DeltaEpsilonVp *= DeltaTime*Sgn(YieldViscoPlastic)/viscosity;
        	rNewStaticData->mEpsilonVp = OldStaticData->mEpsilonVp + DeltaEpsilonVp;
        	// accumulated inelastic strain and damage
        	double DeltaKappaInelastic;
        	DeltaKappaInelastic = this->mDamageDistribution * pow(DeltaEpsilonVp.Norm(),2.);
        	DeltaKappaInelastic +=  (1. - this->mDamageDistribution) * pow(DeltaEpsilonP.Norm(),2.);
        	DeltaKappaInelastic = std::sqrt(DeltaKappaInelastic);
            rNewStaticData->mKappaInelastic = OldStaticData->mKappaInelastic + DeltaKappaInelastic;
            rNewStaticData->mOmegaCompr = 1. - exp(-rNewStaticData->mKappaInelastic/CalculateDuctility());    // additional parameter epsC is neccessary for exp(-mKappaInelastic/epsC);
        }

        // calculate algorithmic tangent
        if (rNewTangent!=0) {
            NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rMatrixMultipl;

            rMatrixMultipl = -((this->DResidualAn(Parameter,Unknown)).fullPivLu().solve(this->DResidualDEpsAn(Unknown))).eval();
            (*rNewTangent) = rMatrixMultipl.block<6,6>(0,0);
//            (*rNewTangent) = this->DStressDEpsNum(ParameterList,Unknown,6); // use for numerical algorithmic tangent
        }

    }

	return Error::SUCCESSFUL;
}

// calculate coefficients of the material matrix
void NuTo::DamageViscoPlasticityEngineeringStress::CalculateCoefficients2DPlainStress(double& C11, double& C12, double& C33) const
{
    double factor = this->mE/(1.0 - (this->mNu * this->mNu));
    C11 = factor;
    C12 = factor * this->mNu;
    C33 = factor * 0.5 * (1.0 - this->mNu);
}

// calculate coefficients of the material matrix
void NuTo::DamageViscoPlasticityEngineeringStress::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE/((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE/(2.*(1.0 + this->mNu));
}

// calculate ductility parameter
double NuTo::DamageViscoPlasticityEngineeringStress::CalculateDuctility()const
{
    return 2.*mFractureEnergy/mCompressiveStrength;
}
///////////////////////////////////////////////////////////////////////////

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get density
//! @return ... density
double NuTo::DamageViscoPlasticityEngineeringStress::GetDensity() const
{
	return this->mRho;
}

//! @brief ... set density
//! @param rRho ... density
void NuTo::DamageViscoPlasticityEngineeringStress::SetDensity(double rRho)
{
    this->CheckDensity(rRho);
    this->mRho = rRho;
    this->SetParametersValid();
}

//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::DamageViscoPlasticityEngineeringStress::GetYoungsModulus() const
{
	return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::DamageViscoPlasticityEngineeringStress::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::DamageViscoPlasticityEngineeringStress::GetPoissonsRatio() const
{
    return mNu;
}


//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::DamageViscoPlasticityEngineeringStress::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}

//! @brief ... get thermal expansion coefficient
//! @return ... thermal expansion coefficient
double NuTo::DamageViscoPlasticityEngineeringStress::GetThermalExpansionCoefficient() const
{
    return mThermalExpansionCoefficient;
}

//! @brief ... set thermal expansion coefficient
//! @param rNu ... thermal expansion coefficient
void NuTo::DamageViscoPlasticityEngineeringStress::SetThermalExpansionCoefficient(double rAlpha)
{
    this->CheckThermalExpansionCoefficient(rAlpha);
    this->mThermalExpansionCoefficient = rAlpha;
    this->SetParametersValid();
}

//! @brief ... get tensile strength
//! @return ... tensile strength
double NuTo::DamageViscoPlasticityEngineeringStress::GetTensileStrength() const
{
    return mTensileStrength;
}

//! @brief ... set tensile strength
//! @param rTensileStrength...  tensile strength
void NuTo::DamageViscoPlasticityEngineeringStress::SetTensileStrength(double rTensileStrength)
{
    this->CheckTensileStrength(rTensileStrength);
    this->mTensileStrength = rTensileStrength;
    this->SetParametersValid();
}

//! @brief ... get compressive strength
//! @return ... compressive strength
double NuTo::DamageViscoPlasticityEngineeringStress::GetCompressiveStrength() const
{
    return mCompressiveStrength;
}

//! @brief ... set compressive strength
//! @param rCompressiveStrength...  compressive strength
void NuTo::DamageViscoPlasticityEngineeringStress::SetCompressiveStrength(double rCompressiveStrength)
{
    this->CheckCompressiveStrength(rCompressiveStrength);
    this->mCompressiveStrength = rCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get biaxial compressive strength
//! @return ... biaxial compressive strength
double NuTo::DamageViscoPlasticityEngineeringStress::GetBiaxialCompressiveStrength() const
{
    return mBiaxialCompressiveStrength;
}

//! @brief ... set biaxial compressive strength
//! @param rBiaxialCompressiveStrength...  biaxial compressive strength
void NuTo::DamageViscoPlasticityEngineeringStress::SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength)
{
    this->CheckBiaxialCompressiveStrength(rBiaxialCompressiveStrength);
    this->mBiaxialCompressiveStrength = rBiaxialCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get viscosity
//! @return ... viscosity
double NuTo::DamageViscoPlasticityEngineeringStress::GetViscosity() const
{
    return mViscosity;
}

//! @brief ... set viscosity
//! @param rViscosity...  viscosity
void NuTo::DamageViscoPlasticityEngineeringStress::SetViscosity(double rViscosity)
{
    this->CheckViscosity(rViscosity);
    this->mViscosity = rViscosity;
    this->SetParametersValid();
}

//! @brief ... get damage distribution (determines the portion of damage via viscoplasticity and plasticity)
//! @return ... damage distribution
double NuTo::DamageViscoPlasticityEngineeringStress::GetDamageDistribution() const
{
    return mDamageDistribution;
}

//! @brief ... set damage distribution (determines the portion of damage via viscoplasticity and plasticity)
//! @param rDamageDistribution... damage distribution
void NuTo::DamageViscoPlasticityEngineeringStress::SetDamageDistribution(double rDamageDistribution)
{
    this->CheckDamageDistribution(rDamageDistribution);
    this->mDamageDistribution = rDamageDistribution;
    this->SetParametersValid();
}

//! @brief ... get viscoplastic yield surface offset with respect to the plastic yield surface
//! @return ... viscoplastic yield surface offset
double NuTo::DamageViscoPlasticityEngineeringStress::GetViscoplasticYieldSurfaceOffset() const
{
    return mViscoplasticYieldSurfaceOffset;
}

//! @brief ... set viscoplastic yield surface offset with respect to the plastic yield surface
//! @param rViscoplasticYieldSurfaceOffset... viscoplastic yield surface offset
void NuTo::DamageViscoPlasticityEngineeringStress::SetViscoplasticYieldSurfaceOffset(double rViscoplasticYieldSurfaceOffset)
{
    this->CheckViscoplasticYieldSurfaceOffset(rViscoplasticYieldSurfaceOffset);
    this->mViscoplasticYieldSurfaceOffset = rViscoplasticYieldSurfaceOffset;
    this->SetParametersValid();
}

//! @brief ... get fracture energy
//! @return ... fracture energy
double NuTo::DamageViscoPlasticityEngineeringStress::GetFractureEnergy() const
{
    return mFractureEnergy;
}

//! @brief ... set fracture energy
//! @param rFractureEnergy... fracture energy
void NuTo::DamageViscoPlasticityEngineeringStress::SetFractureEnergy(double rFractureEnergy)
{
    this->CheckFractureEnergy(rFractureEnergy);
    this->mFractureEnergy = rFractureEnergy;
    this->SetParametersValid();
}

///////////////////////////////////////////////////////////////////////////

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::DamageViscoPlasticityEngineeringStress::GetType() const
{
    return NuTo::Constitutive::LINEAR_ELASTIC;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::DamageViscoPlasticityEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return true;
    case NuTo::Element::PLANE2D10N:
        return true;
    case NuTo::Element::PLANE2D3N:
        return true;
    case NuTo::Element::PLANE2D4N:
        return true;
    case NuTo::Element::PLANE2D4NSPECTRALORDER2:
        return true;
    case NuTo::Element::PLANE2D4NSPECTRALORDER3:
        return true;
    case NuTo::Element::PLANE2D4NSPECTRALORDER4:
        return true;
    case NuTo::Element::PLANE2D6N:
        return true;
    case NuTo::Element::TETRAHEDRON4N:
        return true;
    case NuTo::Element::TETRAHEDRON10N:
        return true;
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if density is positive
//! @param rRho ... density
void NuTo::DamageViscoPlasticityEngineeringStress::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckDensity] The density must be a positive value.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::DamageViscoPlasticityEngineeringStress::CheckYoungsModulus(double rE) const
{
    if (rE < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckYoungsModulus] The Young's modulus must be a non-negative value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::DamageViscoPlasticityEngineeringStress::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check thermal expansion coefficient
//! @param rAlpha ... thermal expansion coefficient
void NuTo::DamageViscoPlasticityEngineeringStress::CheckThermalExpansionCoefficient(double rAlpha) const
{
}

//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ...
void NuTo::DamageViscoPlasticityEngineeringStress::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if compressive strength is positive
//! @param rRadius ... compressive strength
void NuTo::DamageViscoPlasticityEngineeringStress::CheckCompressiveStrength(double rCompressiveStrength) const
{
    if (rCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckCompressiveStrength] The compressive strength must be a positive value.");
    }
}

//! @brief ... check if biaxial compressive strength is positive
//! @param rBiaxialCompressiveStrength ... biaxial compressive strength
void NuTo::DamageViscoPlasticityEngineeringStress::CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const
{
    if (rBiaxialCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be a positive value.");
    }
    if (rBiaxialCompressiveStrength <= mCompressiveStrength)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be higher than the uniaxial compressive strength.");
    }
}

//! @brief ... check whether viscosity is positive
//! @param rViscosity ... viscosity
void NuTo::DamageViscoPlasticityEngineeringStress::CheckViscosity(double rViscosity) const
{
    if (rViscosity <= 0.)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckViscosity] The viscosity must be a positive value.");
    }
}

//! @brief ... check whether the damage distribution ranges between 0 an 1
//! @param rDamageDistribution ... damage distribution
void NuTo::DamageViscoPlasticityEngineeringStress::CheckDamageDistribution(double rDamageDistribution) const
{
    if (rDamageDistribution < 0. || rDamageDistribution > 1.)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckrDamageDistribution] The damage distribution must be between 0 and 1.");
    }
}

//! @brief ... check whether the viscoplastic yield surface offset is negative
//! @param rViscoplasticYieldSurfaceOffset ... viscoplastic yield surface offset
void NuTo::DamageViscoPlasticityEngineeringStress::CheckViscoplasticYieldSurfaceOffset(double rViscoplasticYieldSurfaceOffset) const
{
    if (rViscoplasticYieldSurfaceOffset > 0.)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckrViscoplasticYieldSurfaceOffset] The viscoplastic yield surface offset must be a negative value.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::DamageViscoPlasticityEngineeringStress::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::DamageViscoPlasticityEngineeringStress::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::DamageViscoPlasticityEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus               : " << this->mE << "\n";
    rLogger << "    Poisson's ratio               : " << this->mNu << "\n";
    rLogger << "    Density                       : " << this->mRho << "\n";
    rLogger << "    thermal expansion coefficient : " << this->mThermalExpansionCoefficient << "\n";
    rLogger << "    tensile strength              : " << this->mTensileStrength << "\n";
    rLogger << "    compressive strength          : " << this->mCompressiveStrength << "\n";
    rLogger << "    biaxial compressive strength  : " << this->mBiaxialCompressiveStrength << "\n";
    rLogger << "    viscosity  					  : " << this->mViscosity << "\n";
    rLogger << "    damage distribution  		  : " << this->mDamageDistribution << "\n";
    rLogger << "    viscoplastic yield surface offset	: " << this->mViscoplasticYieldSurfaceOffset << "\n";
    rLogger << "    fracture energy 	 		  : " << this->mFractureEnergy << "\n";
}

// check parameters
void NuTo::DamageViscoPlasticityEngineeringStress::CheckParameters()const
{
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckDensity(this->mRho);
    this->CheckThermalExpansionCoefficient(this->mThermalExpansionCoefficient);
    this->CheckTensileStrength(this->mTensileStrength);
    this->CheckCompressiveStrength(this->mCompressiveStrength);
    this->CheckBiaxialCompressiveStrength(this->mBiaxialCompressiveStrength);
    this->CheckViscosity(this->mViscosity);
    this->CheckDamageDistribution(this->mDamageDistribution);
    this->CheckViscoplasticYieldSurfaceOffset(this->mViscoplasticYieldSurfaceOffset);
    this->CheckFractureEnergy(this->mFractureEnergy);
}

