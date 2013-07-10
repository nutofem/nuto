// $Id: MisesPlasticityEngineeringStress.cpp 112 2009-11-17 16:43:15Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/Logger.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/MisesPlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticity3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticityWithEnergy3D.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/math/FullMatrix.h"

NuTo::MisesPlasticityEngineeringStress::MisesPlasticityEngineeringStress() : ConstitutiveBase()
{
	mE = 0.;
	mNu = 0.;
	mEnergyFlag = false;
	mSigma.resize(1);
	mH.resize(1);
	mThermalExpansionCoefficient = 0;
	SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::MisesPlasticityEngineeringStress::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::MisesPlasticityEngineeringStress::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    rLogger << "start serialize MisesPlasticityEngineeringStress" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
       & BOOST_SERIALIZATION_NVP(mE)
       & BOOST_SERIALIZATION_NVP(mNu)
       & BOOST_SERIALIZATION_NVP(mSigma)
       & BOOST_SERIALIZATION_NVP(mH)
       & BOOST_SERIALIZATION_NVP(mEnergyFlag)
       & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient);
#ifdef DEBUG_SERIALIZATION
    rLogger << "finish serialize MisesPlasticityEngineeringStress" << "\n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::MisesPlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::MisesPlasticityEngineeringStress::Evaluate1D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::Evaluate1D] not implemented for 1D.");
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::MisesPlasticityEngineeringStress::Evaluate2D(ElementBase* rElement, int rIp,
		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::Evaluate2D] not implemented for 2D.");
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::MisesPlasticityEngineeringStress::Evaluate3D(ElementBase* rElement, int rIp,
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

    //Calulate strain
    EngineeringStrain3D engineeringStrain;
	// calculate engineering strain
    auto itDeformationGradient = rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D);
	if(itDeformationGradient==rConstitutiveInput.end())
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::Evaluate3D] deformation gradient 3d needed to evaluate engineering strain3d.");
	const DeformationGradient3D& deformationGradient(itDeformationGradient->second->GetDeformationGradient3D());
	deformationGradient.GetEngineeringStrain(engineeringStrain);

	//substract thermal strain
	EngineeringStrain3D elasticEngineeringStrain(engineeringStrain);
	// if temperature is an input, subtract thermal strains to get elastic strains
	if (section->GetInputConstitutiveIsTemperature())
	{
		std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
		if (itInput==rConstitutiveInput.end())
			throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::Evaluate3D] temperature needed to evaluate thermal engineering strain3d.");
		double temperature(itInput->second->GetTemperature());
		double deltaStrain(mThermalExpansionCoefficient * temperature);
		EngineeringStrain3D elasticEngineeringStrain;
		elasticEngineeringStrain[0] -= deltaStrain;
		elasticEngineeringStrain[1] -= deltaStrain;
		elasticEngineeringStrain[2] -= deltaStrain;
	}

	EngineeringStress3D engineeringStress;
	ConstitutiveTangentLocal<6,6>* tangent(0);
	ConstitutiveStaticDataMisesPlasticity3D newStaticData;

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

	//check, if an nonlinear iteration has to be performed, in this simple case, just calculate the linear elastic coefficients
    if (performReturnMapping)
    {
        // perform return mapping
    	NuTo::Error::eError errorReturnMapping = ReturnMapping3D(rElement, rIp, elasticEngineeringStrain, &engineeringStress, tangent, &newStaticData ,rElement->GetStructure()->GetLogger());
    	if (errorReturnMapping!=Error::SUCCESSFUL)
    		return errorReturnMapping;
    }

    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
    		itOutput != rConstitutiveOutput.end(); itOutput++)
    {
    	switch(itOutput->first)
    	{
    	case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
    	{
			itOutput->second->GetEngineeringStress3D() = engineeringStress;
   		break;
    	}
    	case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D:
    	{
    		//tangent has already been calculated
			itOutput->second->AsConstitutiveTangentLocal_6x6().SetSymmetry(true);

    		break;
    	}
    	case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
    		itOutput->second->GetEngineeringStrain3D() = engineeringStrain;
    		break;
    	case NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D:
    	{
    		EngineeringStrain3D& engineeringPlasticStrain(itOutput->second->GetEngineeringStrain3D());
    		engineeringPlasticStrain[0] = newStaticData.mEpsilonP[0] ;
    		engineeringPlasticStrain[1] = newStaticData.mEpsilonP[1] ;
    		engineeringPlasticStrain[2] = newStaticData.mEpsilonP[2] ;
    		engineeringPlasticStrain[3] = newStaticData.mEpsilonP[3] ;
    		engineeringPlasticStrain[4] = newStaticData.mEpsilonP[4] ;
    		engineeringPlasticStrain[5] = newStaticData.mEpsilonP[5] ;
    		break;
    	}
    	case NuTo::Constitutive::Output::DAMAGE:
    	{
    		itOutput->second->GetDamage().SetDamage(0.);
   		break;
    	}
    	case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
    	break;
    	case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
    	{
    	    *(rElement->GetStaticData(rIp)->AsConstitutiveStaticDataMisesPlasticity3D()) = newStaticData;
    	}
		break;
    	default:
    		throw MechanicsException(std::string("[NuTo::MisesPlasticityEngineeringStress::Evaluate3D] output object)") +
    				NuTo::Constitutive::OutputToString(itOutput->first) +
    				std::string(" culd not be calculated, check the allocated material law and the section behavior."));
    	}
    }

    //update history variables but for linear elastic, there is nothing to do

	return Error::SUCCESSFUL;
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::MisesPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
		const ElementBase* rElement) const
{
	throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D] To be implemented.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::MisesPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
		const ElementBase* rElement) const
{
	throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D] To be implemented.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::MisesPlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
		const ElementBase* rElement) const
{
	if (mEnergyFlag)
	{
		return new NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D();
	}
	else
	{
		return new NuTo::ConstitutiveStaticDataMisesPlasticity3D();
	}
}

//! @brief ... performs the return mapping procedure in 3D
//! @param rElement ... structure
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient point
//! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
//! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
//! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
#define sqrt_2div3 0.81649658
#define tolerance 1e-8
NuTo::Error::eError NuTo::MisesPlasticityEngineeringStress::ReturnMapping3D(const ElementBase* rElement,int rIp,
		const EngineeringStrain3D& rEngineeringStrain,
		EngineeringStress3D* rNewStress,
		ConstitutiveTangentLocal<6,6>* rNewTangent,
		ConstitutiveStaticDataMisesPlasticity3D* rNewStaticData,
		Logger& rLogger)const
{
    double sigma_trial[6],
    xi_trial[6],
    norm_dev,
    sigma_y,
    factor,
    factor2,
    yield_condition,
    epsilon_p_eq2,
    d_sigma,
    d_H,
    H,
    H2,
    mu,
    bulk_modulus,
    delta_gamma=0.,
    g,
    dg,
    df_dsigma[6],
    trace_epsilon,
    trace_epsilon_div_3;

    //here use the modified parameters to enable random fields
    double modE = GetRanfieldFactorYoungsModulus(rElement,rIp)*mE;
    CheckYoungsModulus(modE);

    double modNu = GetRanfieldFactorPoissonsRatio(rElement,rIp)*mNu;
    CheckPoissonsRatio(modNu);

    double ranfield_factor_yield_strength(GetRanfieldFactorYieldStrength(rElement,rIp));
    double ranfield_factor_hardening_modulus(GetRanfieldFactorHardeningModulus(rElement,rIp));

    mu    = modE/(2.*(1.+modNu));
    bulk_modulus = modE/(3.-6.*modNu);

    // set strain data ptr
    const double *total_strain(&(rEngineeringStrain[0]));

    //get old static data
    const ConstitutiveStaticDataMisesPlasticity3D* rOldStaticData = rElement->GetStaticData(rIp)->AsConstitutiveStaticDataMisesPlasticity3D();
    const double *plastic_strain(rOldStaticData->mEpsilonP);
    const double *back_stress(rOldStaticData->mSigmaB);

    trace_epsilon = total_strain[0] + total_strain[1] + total_strain[2];

    trace_epsilon_div_3 = trace_epsilon/3.;

    //trial stress
    sigma_trial[0] = (total_strain[0]-trace_epsilon_div_3-plastic_strain[0])*2.*mu;
    sigma_trial[1] = (total_strain[1]-trace_epsilon_div_3-plastic_strain[1])*2.*mu;
    sigma_trial[2] = (total_strain[2]-trace_epsilon_div_3-plastic_strain[2])*2.*mu;
    sigma_trial[3] = (total_strain[3]                    -plastic_strain[3])*mu; // in total strain, gamma is stored
    sigma_trial[4] = (total_strain[4]                    -plastic_strain[4])*mu; // in total strain, gamma is stored
    sigma_trial[5] = (total_strain[5]                    -plastic_strain[5])*mu; // in total strain, gamma is stored

    //subtract backstress
    xi_trial[0] = sigma_trial[0]-back_stress[0];
    xi_trial[1] = sigma_trial[1]-back_stress[1];
    xi_trial[2] = sigma_trial[2]-back_stress[2];
    xi_trial[3] = sigma_trial[3]-back_stress[3];
    xi_trial[4] = sigma_trial[4]-back_stress[4];
    xi_trial[5] = sigma_trial[5]-back_stress[5];

    //printf("total_strain %g %g %g %g %g %g\n",total_strain[0],total_strain[1],total_strain[2],total_strain[3],total_strain[4],total_strain[5]);
    //printf("Xi trial %g %g %g %g %g %g\n",xi_trial[0],xi_trial[1],xi_trial[2],xi_trial[3],xi_trial[4],xi_trial[5]);

    // norm of deviator
    norm_dev = sqrt(xi_trial[0]*xi_trial[0]+xi_trial[1]*xi_trial[1]+xi_trial[2]*xi_trial[2]+
                    2.*(xi_trial[3]*xi_trial[3]+xi_trial[4]*xi_trial[4]+xi_trial[5]*xi_trial[5]));

    //determine radius of yield function
    sigma_y = ranfield_factor_yield_strength*GetYieldStrength(rOldStaticData->mEpsilonPEq,d_sigma);
    d_sigma*= ranfield_factor_yield_strength;
    yield_condition = norm_dev - sqrt_2div3 * sigma_y;

    if (yield_condition<-tolerance*sigma_y)
    {
        // elastic regime
    	factor = bulk_modulus*trace_epsilon;
        if (rNewStress!=0)
        {
        	(*rNewStress)[0] = factor+sigma_trial[0];
        	(*rNewStress)[1] = factor+sigma_trial[1];
        	(*rNewStress)[2] = factor+sigma_trial[2];
        	(*rNewStress)[3] = 	   sigma_trial[3];
        	(*rNewStress)[4] = 	   sigma_trial[4];
        	(*rNewStress)[5] = 	   sigma_trial[5];
        }
        if (rNewTangent!=0)
        {
            factor = modE/(1.+modNu)/(1.-2.*modNu);

            (*rNewTangent)(0,0) = (1.-modNu)*factor;
            (*rNewTangent)(1,0) = modNu*factor;
            (*rNewTangent)(2,0) = (*rNewTangent)(1,0);
            (*rNewTangent)(3,0) = 0.;
            (*rNewTangent)(4,0) = 0.;
            (*rNewTangent)(5,0) = 0.;

            (*rNewTangent)(0,1) = (*rNewTangent)(1,0);
            (*rNewTangent)(1,1) = (*rNewTangent)(0,0);
            (*rNewTangent)(2,1) = (*rNewTangent)(1,0);
            (*rNewTangent)(3,1) = 0.;
            (*rNewTangent)(4,1) = 0.;
            (*rNewTangent)(5,1) = 0.;

            (*rNewTangent)(0,2) = (*rNewTangent)(1,0);
            (*rNewTangent)(1,2) = (*rNewTangent)(1,0);
            (*rNewTangent)(2,2) = (*rNewTangent)(0,0);
            (*rNewTangent)(3,2) = 0.;
            (*rNewTangent)(4,2) = 0.;
            (*rNewTangent)(5,2) = 0.;

            (*rNewTangent)(0,3) = 0.;
            (*rNewTangent)(1,3) = 0.;
            (*rNewTangent)(2,3) = 0.;
            (*rNewTangent)(3,3) = (0.5-modNu)*factor;
            (*rNewTangent)(4,3) = 0.;
            (*rNewTangent)(5,3) = 0.;

            (*rNewTangent)(0,4) = 0.;
            (*rNewTangent)(1,4) = 0.;
            (*rNewTangent)(2,4) = 0.;
            (*rNewTangent)(3,4) = 0.;
            (*rNewTangent)(4,4) = (*rNewTangent)(3,3);
            (*rNewTangent)(5,4) = 0.;

            (*rNewTangent)(0,5) = 0.;
            (*rNewTangent)(1,5) = 0.;
            (*rNewTangent)(2,5) = 0.;
            (*rNewTangent)(3,5) = 0.;
            (*rNewTangent)(4,5) = 0.;
            (*rNewTangent)(5,5) = (*rNewTangent)(3,3);
        }

        // static data is unchanged
        return Error::SUCCESSFUL;
    }

    //plastic loading
    H  = ranfield_factor_hardening_modulus * GetHardeningModulus(rOldStaticData->mEpsilonPEq,d_H);
    d_H *= ranfield_factor_hardening_modulus;
    epsilon_p_eq2 = rOldStaticData->mEpsilonPEq;
    H2 = H;

    int i=0;
    for (;i<100;i++)
    {
        g  = yield_condition - (2.*mu*delta_gamma + sqrt_2div3 * (H2-H));
        if (fabs(g)<tolerance*sigma_y)
        {
            break;
        }
        dg = -2.* mu * (1.+(d_H+d_sigma)/(3.*mu));
        delta_gamma -=g/dg;
        epsilon_p_eq2 = rOldStaticData->mEpsilonPEq + sqrt_2div3 * delta_gamma;

        H2  = ranfield_factor_hardening_modulus * GetHardeningModulus(epsilon_p_eq2,d_H);
        d_H *= ranfield_factor_hardening_modulus;

        sigma_y = ranfield_factor_yield_strength*GetYieldStrength(epsilon_p_eq2,d_sigma);
        d_sigma*= ranfield_factor_yield_strength;

        yield_condition = norm_dev - sqrt_2div3 * sigma_y;
    }

    if (i==100)
    {
        rLogger << "yield condition " << yield_condition << " delta_gamma " << delta_gamma << "\n";
        rLogger << "epsilon_p_eq " << rOldStaticData->mEpsilonPEq;
        rLogger << "total strain " << total_strain[0] << " " << total_strain[1] << " " << total_strain[2] << " " << total_strain[3] << " " << total_strain[4] << " " << total_strain[5] << " " <<"\n";
        rLogger << "plastic strain " << plastic_strain[0] << " " << plastic_strain[1] << " " << plastic_strain[2] << " " << plastic_strain[3] << " " << plastic_strain[4] << " " << plastic_strain[5] << " " <<"\n";
        rLogger << "back stress" << back_stress[0] << " " << back_stress[1] << " " << back_stress[2] << " " << back_stress[3] << " " << back_stress[4] << " " << back_stress[5] << " " <<"\n";
        rLogger << "[NuTo::MisesPlasticityEngineeringStress::ReturnMapping3D] No convergence after 100 steps, check the source code." << "\n";

        return Error::NO_CONVERGENCE;

    }

    /* derivative of yield surface */
    df_dsigma[0]  = xi_trial[0]/norm_dev;
    df_dsigma[1]  = xi_trial[1]/norm_dev;
    df_dsigma[2]  = xi_trial[2]/norm_dev;
    df_dsigma[3]  = xi_trial[3]/norm_dev;
    df_dsigma[4]  = xi_trial[4]/norm_dev;
    df_dsigma[5]  = xi_trial[5]/norm_dev;

    //update static data
    if (rNewStaticData!=0)
    {
    	//update equivalent plastic strain
    	rNewStaticData->mEpsilonPEq = rOldStaticData->mEpsilonPEq + sqrt_2div3 * delta_gamma;

        //update backstress
        factor = sqrt_2div3 * (H2-H);
        rNewStaticData->mSigmaB[0] = rOldStaticData->mSigmaB[0] + factor*df_dsigma[0];
        rNewStaticData->mSigmaB[1] = rOldStaticData->mSigmaB[1] + factor*df_dsigma[1];
        rNewStaticData->mSigmaB[2] = rOldStaticData->mSigmaB[2] + factor*df_dsigma[2];
        rNewStaticData->mSigmaB[3] = rOldStaticData->mSigmaB[3] + factor*df_dsigma[3];
        rNewStaticData->mSigmaB[4] = rOldStaticData->mSigmaB[4] + factor*df_dsigma[4];
        rNewStaticData->mSigmaB[5] = rOldStaticData->mSigmaB[5] + factor*df_dsigma[5];

        //update plastic_strain
        rNewStaticData->mEpsilonP[0] = rOldStaticData->mEpsilonP[0] + delta_gamma*df_dsigma[0];
        rNewStaticData->mEpsilonP[1] = rOldStaticData->mEpsilonP[1] + delta_gamma*df_dsigma[1];
        rNewStaticData->mEpsilonP[2] = rOldStaticData->mEpsilonP[2] + delta_gamma*df_dsigma[2];
        rNewStaticData->mEpsilonP[3] = rOldStaticData->mEpsilonP[3] + 2.*delta_gamma*df_dsigma[3];  /* gamma */
        rNewStaticData->mEpsilonP[4] = rOldStaticData->mEpsilonP[4] + 2.*delta_gamma*df_dsigma[4];  /* gamma */
        rNewStaticData->mEpsilonP[5] = rOldStaticData->mEpsilonP[5] + 2.*delta_gamma*df_dsigma[5];  /* gamma */
    }

    //update stress
    if (rNewStress!=0)
    {
		factor  = 2.*mu*delta_gamma;
		factor2 = bulk_modulus*trace_epsilon;
		(*rNewStress)[0] = factor2+sigma_trial[0]-factor*df_dsigma[0];
		(*rNewStress)[1] = factor2+sigma_trial[1]-factor*df_dsigma[1];
		(*rNewStress)[2] = factor2+sigma_trial[2]-factor*df_dsigma[2];
		(*rNewStress)[3] =         sigma_trial[3]-factor*df_dsigma[3];
		(*rNewStress)[4] =         sigma_trial[4]-factor*df_dsigma[4];
		(*rNewStress)[5] =         sigma_trial[5]-factor*df_dsigma[5];
    }

    //update stiffness
    if (rNewTangent!=0)
    {
        double theta     = 1.-2.*mu*delta_gamma/norm_dev;
        double theta_bar = 1./(1.+(d_sigma+d_H)/(3.*mu))-(1.-theta);
        factor    = 2.*mu*theta;
        double factor_div3  = factor/3.;
        double factor_mul2div3  = 2.*factor_div3;
        double factor2   = -2.*mu*theta_bar;
        double factor3   = bulk_modulus;

        (*rNewTangent)(0,0) =   (factor3 + factor_mul2div3 + factor2*df_dsigma[0]*df_dsigma[0]);
        (*rNewTangent)(1,0) =   (factor3 - factor_div3     + factor2*df_dsigma[0]*df_dsigma[1]);
        (*rNewTangent)(2,0) =   (factor3 - factor_div3     + factor2*df_dsigma[0]*df_dsigma[2]);
        (*rNewTangent)(3,0) =   (							 factor2*df_dsigma[0]*df_dsigma[3]);
        (*rNewTangent)(4,0) =   ( 						 factor2*df_dsigma[0]*df_dsigma[4]);
        (*rNewTangent)(5,0) =   ( 						 factor2*df_dsigma[0]*df_dsigma[5]);
        (*rNewTangent)(0,1) =   (*rNewTangent)(1,0);
        (*rNewTangent)(1,1) =   (factor3 + factor_mul2div3 + factor2*df_dsigma[1]*df_dsigma[1]);
        (*rNewTangent)(2,1) =   (factor3 - factor_div3     + factor2*df_dsigma[1]*df_dsigma[2]);
        (*rNewTangent)(3,1) =   (						 factor2*df_dsigma[1]*df_dsigma[3]);
        (*rNewTangent)(4,1) =   (						 factor2*df_dsigma[1]*df_dsigma[4]);
        (*rNewTangent)(5,1) =   (						 factor2*df_dsigma[1]*df_dsigma[5]);
        (*rNewTangent)(0,2) =   (*rNewTangent)(2,0);
        (*rNewTangent)(1,2) =   (*rNewTangent)(2,1);
        (*rNewTangent)(2,2) =   (factor3 + factor_mul2div3 + factor2*df_dsigma[2]*df_dsigma[2]);
        (*rNewTangent)(3,2) =   (						 factor2*df_dsigma[2]*df_dsigma[3]);
        (*rNewTangent)(4,2) =   (						 factor2*df_dsigma[2]*df_dsigma[4]);
        (*rNewTangent)(5,2) =   (						 factor2*df_dsigma[2]*df_dsigma[5]);
        (*rNewTangent)(0,3) =   (*rNewTangent)(3,0);
        (*rNewTangent)(1,3) =   (*rNewTangent)(3,1);
        (*rNewTangent)(2,3) =   (*rNewTangent)(3,2);
        (*rNewTangent)(3,3) =   (	 0.5*factor 		 +factor2*df_dsigma[3]*df_dsigma[3]);
        (*rNewTangent)(4,3) =   (					      factor2*df_dsigma[3]*df_dsigma[4]);
        (*rNewTangent)(5,3) =   (					      factor2*df_dsigma[3]*df_dsigma[5]);
        (*rNewTangent)(0,4) =   (*rNewTangent)(4,0);
        (*rNewTangent)(1,4) =   (*rNewTangent)(4,1);
        (*rNewTangent)(2,4) =   (*rNewTangent)(4,2);
        (*rNewTangent)(3,4) =   (*rNewTangent)(4,3);
        (*rNewTangent)(4,4) =   (	 0.5*factor 		 +factor2*df_dsigma[4]*df_dsigma[4]);
        (*rNewTangent)(5,4) =   (					      factor2*df_dsigma[4]*df_dsigma[5]);
        (*rNewTangent)(0,5) =   (*rNewTangent)(5,0);
        (*rNewTangent)(1,5) =   (*rNewTangent)(5,1);
        (*rNewTangent)(2,5) =   (*rNewTangent)(5,2);
        (*rNewTangent)(3,5) =   (*rNewTangent)(5,3);
        (*rNewTangent)(4,5) =   (*rNewTangent)(5,4);
        (*rNewTangent)(5,5) =   (	 0.5*factor 		 +factor2*df_dsigma[5]*df_dsigma[5]);
        rNewTangent->SetSymmetry(true);
    }

    return Error::SUCCESSFUL;
}


//! @brief ... calculates for a given equivalent plastic strain the radius of the yield surface
//! @param ... rEpsilonPEq equivalent plastic strain
//! @param ... rDSigmaDEpsilonP derivative of the yield strength with respect to the plastic strains (return value)
//! @return ... yield strength (radius of the yield surface)
double NuTo::MisesPlasticityEngineeringStress::GetYieldStrength(double rEpsilonPEq, double& rDSigmaDEpsilonP)const
{
    assert(mSigma.size()>0);
    std::vector<std::pair<double, double> >::const_iterator it(mSigma.begin());
    while (rEpsilonPEq>=it->first)
    {
        it++;
        if (it==mSigma.end())
        {
            // the maximum is reached, afterwards the yield strength remains constant
        	rDSigmaDEpsilonP = 0.;
            return (it-1)->second;
        }
    }

  	rDSigmaDEpsilonP = (it->second-(it-1)->second)/(it->first-(it-1)->first);
   	return rDSigmaDEpsilonP*(rEpsilonPEq-(it-1)->first)+(it-1)->second;
}

//! @brief ... calculates for a given equivalent plastic strain the hardening modulus
//! @param ... rEpsilonPEq equivalent plastic strain
//! @param ... rDSigmaDEpsilonP derivative of the hardening modulus with respect to the plastic strains (return value)
//! @return ... hardening modulus
double NuTo::MisesPlasticityEngineeringStress::GetHardeningModulus(double rEpsilonPEq, double& rDHDEpsilonP)const
{
    assert(mH.size()>0);
	std::vector<std::pair<double, double> >::const_iterator it(mH.begin());
    double H(0);
    if (mH.size()==1)
    {
        rDHDEpsilonP =it->second;
    	H=rEpsilonPEq*rDHDEpsilonP;
        return H;
    }
    else
    {
		do
		{
			it++;
			if (it==mH.end())
			{
				// the maximum is reached, afterwards use a constant slope
				it--;
				rDHDEpsilonP = it->second;
				H+=(rEpsilonPEq-it->first)*(it->second);
				return H;
			}
			if (rEpsilonPEq>it->first)
			{
				H+=((it)->first-(it-1)->first)*(it-1)->second;
			}
			else
			{
				it--;
				rDHDEpsilonP = it->second;
				H+=(rEpsilonPEq-it->first)*(it->second);
				return H;
			}
		}
		while(true);
    }
}

///////////////////////////////////////////////////////////////////////////

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::MisesPlasticityEngineeringStress::GetYoungsModulus() const
{
	return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::MisesPlasticityEngineeringStress::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::MisesPlasticityEngineeringStress::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::MisesPlasticityEngineeringStress::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}
//! @brief ... get initial yield strength
//! @return ... yield strength
double NuTo::MisesPlasticityEngineeringStress::GetInitialYieldStrength() const
{
	if(mSigma.size()==0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::GetInitialYieldStrength] Size of yield strength vector is zero.");
	return mSigma[0].second;
}

//! @brief ... set initial yield strength
//! @param rSigma ...  yield strength
void NuTo::MisesPlasticityEngineeringStress::SetInitialYieldStrength(double rSigma)
{
	if(mSigma.size()==0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::SetInitialYieldStrength] Size of yield strength vector is zero.");
	if (rSigma<=0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::SetInitialYieldStrength] Initial yield strength has to be positive.");
	mSigma[0].second = rSigma;
    this->SetParametersValid();
}

//! @brief ... get yield strength for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding yield strength
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::MisesPlasticityEngineeringStress::GetYieldStrength() const
{
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> returnMatrix(mSigma.size(),2);
	for (unsigned int count=0; count<mSigma.size(); count++)
	{
		returnMatrix(count,0) = mSigma[count].first;
		returnMatrix(count,1) = mSigma[count].second;
	}
	return returnMatrix;
}

//! @brief ... add yield strength
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  yield strength
void NuTo::MisesPlasticityEngineeringStress::AddYieldStrength(double rEpsilon, double rSigma)
{
	if (rEpsilon<=0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::AddYieldStrength] Equivalente strain has to be positive.");
	if (rSigma<=0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::AddYieldStrength] Yield strength has to be positive.");
	std::vector<std::pair<double, double> >::iterator it;
	for (it = mSigma.begin(); it!=mSigma.end(); it++)
	{
		if (it->first>rEpsilon)
		{
			if (it!=mSigma.begin())
			{
				if ((it-1)->second>(it->second))
					throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::AddYieldStrength] The yield strength can only increase for increasing epsilon equivalent.");
			}
			break;
		}
	}
	mSigma.insert(it,1,std::pair<double,double>(rEpsilon,rSigma));
    this->SetParametersValid();
}

//! @brief ... get initial hardening modulus
//! @return ... hardening modulus
double NuTo::MisesPlasticityEngineeringStress::GetInitialHardeningModulus() const
{
	if(mH.size()==0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::SetInitialHardeningModulus] Size of hardening modulus vector is zero.");
	return mH[0].second;
}

//! @brief ... set initial hardening modulus
//! @param rH ...  hardening modulus
void NuTo::MisesPlasticityEngineeringStress::SetInitialHardeningModulus(double rH)
{
	if(mH.size()==0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::SetInitialHardeningModulus] Size of hardening modulus vector is zero.");
	if (rH<0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::SetInitialHardeningModulus] Initial hardening modulus must not be negative.");
	mH[0].second = rH;
    this->SetParametersValid();
}

//! @brief ... get hardening modulus for multilinear response
//! @return ... first column: equivalent plastic strain
//! @return ... second column: corresponding hardening modulus
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::MisesPlasticityEngineeringStress::GetHardeningModulus() const
{
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> returnMatrix(mH.size(),2);
	for (unsigned int count=0; count<mH.size(); count++)
	{
		returnMatrix(count,0) = mH[count].first;
		returnMatrix(count,1) = mH[count].second;
	}
	return returnMatrix;
}

//! @brief ... add hardening modulus
//! @param rEpsilon ...  equivalent plastic strain
//! @param rSigma ...  hardening modulus
void NuTo::MisesPlasticityEngineeringStress::AddHardeningModulus(double rEpsilon, double rH)
{
	if (rEpsilon<=0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::AddHardeningModulus] Equivalente strain has to be positive.");
	if (rH<0)
		throw MechanicsException("[NuTo::MisesPlasticityEngineeringStress::AddHardeningModulus] Hardening modul must not be negative.");
	std::vector<std::pair<double, double> >::iterator it;
	for (it = mH.begin(); it!=mH.end(); it++)
	{
		if (it->first>rEpsilon)
		{
			break;
		}
	}
	mH.insert(it,1,std::pair<double,double>(rEpsilon,rH));
    this->SetParametersValid();
}

//! @brief ... get thermal expansion coefficient
//! @return ... thermal expansion coefficient
double NuTo::MisesPlasticityEngineeringStress::GetThermalExpansionCoefficient() const
{
    return mThermalExpansionCoefficient;
}

//! @brief ... set thermal expansion coefficient
//! @param rNu ... thermal expansion coefficient
void NuTo::MisesPlasticityEngineeringStress::SetThermalExpansionCoefficient(double rAlpha)
{
    this->CheckThermalExpansionCoefficient(rAlpha);
    this->mThermalExpansionCoefficient = rAlpha;
    this->SetParametersValid();
}

///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::MisesPlasticityEngineeringStress::GetType() const
{
    return NuTo::Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::MisesPlasticityEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::TETRAHEDRON4N:
        return true;
    case NuTo::Element::TETRAHEDRON10N:
        return true;
    case NuTo::Element::BRICK8N:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::MisesPlasticityEngineeringStress::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::MisesPlasticityEngineeringStress::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check yield strength is positive
//! @param rSigma ... yield strength
void NuTo::MisesPlasticityEngineeringStress::CheckYieldStrength(std::vector<std::pair<double, double> > rSigma) const
{
	if (rSigma.size()==0)
    {
        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckYieldStrength] At least an initial yield strength is required.");
    }
	if (rSigma[0].second <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckYieldStrength] The initial yield strength must be a positive value.");
    }
	rSigma[0].first = 0.;

	for (unsigned int count=1; count<rSigma.size(); count++)
	{
		if (rSigma[count-1].first>=rSigma[count].first)
	        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckYieldStrength] For multilinear plasticity, the epsilon should always increase.");
		if (rSigma[count-1].second>rSigma[count].second)
	        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckYieldStrength] For multilinear plasticity, the yield strength should always increase.");
	}
}

//! @brief ... check hardening modulus
//! @param rH ... hardening modulus
void NuTo::MisesPlasticityEngineeringStress::CheckHardeningModulus(std::vector<std::pair<double, double> > rH) const
{
	if (rH.size()==0)
    {
        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckHardeningModulus] At least an initial hardening modulus is required.");
    }
	if (rH[0].second < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckHardeningModulus] The initial hardening modulus must not be a negative value.");
    }
	rH[0].first = 0.;

	for (unsigned int count=1; count<rH.size(); count++)
	{
		if (rH[count-1].first>=rH[count].first)
	        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckHardeningModulus] For multilinear plasticity, the epsilon should always increase.");
		if (rH[count].second<0)
	        throw NuTo::MechanicsException("[NuTo::MisesPlasticityEngineeringStress::CheckHardeningModulus] For multilinear plasticity, the hardening modulus should always be positive.");
	}
}

//! @brief ... check thermal expansion coefficient
//! @param rAlpha ... thermal expansion coefficient
void NuTo::MisesPlasticityEngineeringStress::CheckThermalExpansionCoefficient(double rAlpha) const
{
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::MisesPlasticityEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus: " << this->mE << "\n";
    rLogger << "    Poisson's ratio: " << this->mNu << "\n";
	rLogger << "    multilinear yield strength: (interval,epsilon,sigma)" << "\n";
    for (unsigned int count=0; count<mSigma.size(); count++)
    	rLogger << "       " << count<< " : " << this->mSigma[count].first << "    " << this->mSigma[count].second << "\n";
	rLogger << "    multilinear hardening modulus: (interval,epsilon,H')" << "\n";
    for (unsigned int count=0; count<mH.size(); count++)
    	rLogger << "       " << count<< " : " << this->mH[count].first << "    " << this->mH[count].second << "\n";
    rLogger << "    mThermalExpansionCoefficient: " << this->mThermalExpansionCoefficient << "\n";
}

// check parameters
void NuTo::MisesPlasticityEngineeringStress::CheckParameters()const
{
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckYieldStrength(this->mSigma);
    this->CheckHardeningModulus(this->mH);
    this->CheckThermalExpansionCoefficient(this->mThermalExpansionCoefficient);
}
