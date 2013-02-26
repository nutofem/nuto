// $Id: GradientDamagePlasticityEngineeringStress.cpp 612 2012-08-13 07:31:23Z unger3 $
// GradientDamagePlasticityEngineeringStress.cpp
// created Apr 26, 2010 by Joerg F. Unger
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <eigen3/Eigen/LU>

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamagePlasticity1D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamagePlasticityEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamage.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#define sqrt3 1.732050808
#define MAX_OMEGA 0.999
#define ENABLE_DEBUG

NuTo::GradientDamagePlasticityEngineeringStress::GradientDamagePlasticityEngineeringStress() : ConstitutiveBase()
{
    mRho = 0.;
    mE = 0.;
    mNu = 0.;
    mNonlocalRadius = 1.;
    mTensileStrength = 0.;
    mCompressiveStrength = 0.;
    mBiaxialCompressiveStrength = 0.;
    mFractureEnergy = 0.;
    mEpsilonF = 1.;
    mYieldSurface = Constitutive::RANKINE_ROUNDED;
    mDamage = true;
    mThermalExpansionCoefficient = 0.;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::GradientDamagePlasticityEngineeringStress::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
       rLogger << "start serialize GradientDamagePlasticityEngineeringStress" << "\n";
#endif
       ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
          & BOOST_SERIALIZATION_NVP(mRho)
          & BOOST_SERIALIZATION_NVP(mE)
          & BOOST_SERIALIZATION_NVP(mNu)
          & BOOST_SERIALIZATION_NVP(mNonlocalRadius)
          & BOOST_SERIALIZATION_NVP(mTensileStrength)
          & BOOST_SERIALIZATION_NVP(mCompressiveStrength)
          & BOOST_SERIALIZATION_NVP(mBiaxialCompressiveStrength)
          & BOOST_SERIALIZATION_NVP(mFractureEnergy)
          & BOOST_SERIALIZATION_NVP(mYieldSurface)
          & BOOST_SERIALIZATION_NVP(mDamage)
          & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient);
#ifdef DEBUG_SERIALIZATION
       rLogger << "finish serialize GradientDamagePlasticityEngineeringStress" << "\n";
#endif
    }
    BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::GradientDamagePlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamagePlasticityEngineeringStress::Evaluate1D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    	// get section information determining which input on the constitutive level should be used
        const SectionBase* section(rElement->GetSection());

        // check if parameters are valid
        if (this->mParametersValid == false)
        {
               //throw an exception giving information related to the wrong parameter
            CheckParameters();
        }

        if (section->GetType()!=Section::TRUSS)
        	throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate] only truss sections are implemented.");

        EngineeringStrain1D engineeringStrain1D;
        // calculate engineering strain
        if(rConstitutiveInput.find(NuTo::Constitutive::eInput::DEFORMATION_GRADIENT_1D)==rConstitutiveInput.end())
            throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate] deformation gradient 1d needed to evaluate engineering strain1d.");
        const DeformationGradient1D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::eInput::DEFORMATION_GRADIENT_1D)->second->GetDeformationGradient1D());
        deformationGradient.GetEngineeringStrain(engineeringStrain1D);

        //Get previous ip_data
        ConstitutiveStaticDataGradientDamagePlasticity1D *oldStaticData = (rElement->GetStaticData(rIp))->AsGradientDamagePlasticity1D();

        //total trial strain in 3D to use the full 3D return mapping
		double prevEpsilonP[6];
		prevEpsilonP[0] = oldStaticData->mEpsilonP[0];
		prevEpsilonP[1] = oldStaticData->mEpsilonP[1];
		prevEpsilonP[2] = oldStaticData->mEpsilonP[1];
		prevEpsilonP[3] = 0.;
		prevEpsilonP[4] = 0.;
		prevEpsilonP[5] = 0.;

        EngineeringStrain3D engineeringStrain3D;
        engineeringStrain3D.mEngineeringStrain[0] = engineeringStrain1D.mEngineeringStrain;
        engineeringStrain3D.mEngineeringStrain[1] = oldStaticData->mEpsilonTotRadial + mNu*(engineeringStrain1D.mEngineeringStrain-oldStaticData->mPrevStrain.mEngineeringStrain); //just use the previous total radial strain plus poisson*delta strain in axial direction
        engineeringStrain3D.mEngineeringStrain[2] = engineeringStrain3D.mEngineeringStrain[1];
        engineeringStrain3D.mEngineeringStrain[3] =  0.;
        engineeringStrain3D.mEngineeringStrain[4] =  0.;
        engineeringStrain3D.mEngineeringStrain[5] =  0.;

        std::cout << "[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate1D] this is just for test purposes." << std::endl;
        engineeringStrain3D.mEngineeringStrain[0] = 1.5*mCompressiveStrength/mE;
        engineeringStrain3D.mEngineeringStrain[1] = 0.0001; //just use the previous total radial strain plus poisson*delta strain in axial direction
        engineeringStrain3D.mEngineeringStrain[2] = 0.0001;

        EngineeringStrain3D prevEngineeringStrain3D;
        prevEngineeringStrain3D.mEngineeringStrain[0] = oldStaticData->mPrevStrain.mEngineeringStrain;
        prevEngineeringStrain3D.mEngineeringStrain[1] = oldStaticData->mEpsilonTotRadial;
        prevEngineeringStrain3D.mEngineeringStrain[2] = oldStaticData->mEpsilonTotRadial;
        prevEngineeringStrain3D.mEngineeringStrain[3] =  0.;
        prevEngineeringStrain3D.mEngineeringStrain[4] =  0.;
        prevEngineeringStrain3D.mEngineeringStrain[5] =  0.;

        // subtract thermal strain
        if (section->GetInputConstitutiveIsTemperature())
        {
            std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::eInput::TEMPERATURE));
            if (itInput==rConstitutiveInput.end())
                throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate2D] temperature needed to evaluate thermal engineering strain2d.");
            double temperature(itInput->second->GetTemperature());
            double deltaStrain(mThermalExpansionCoefficient * temperature);
            engineeringStrain3D.mEngineeringStrain[0] -= deltaStrain;
            engineeringStrain3D.mEngineeringStrain[1] -= deltaStrain;
            engineeringStrain3D.mEngineeringStrain[2] -= deltaStrain;
            throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate1D] add temperature history.");
/*            double prevTemperature(oldStaticData->mPrevTemperature);
            double deltaStrain(mThermalExpansionCoefficient * prevTemperature);
            prevEngineeringStrain3D.mEngineeringStrain[0] -= deltaStrain;
            prevEngineeringStrain3D.mEngineeringStrain[1] -= deltaStrain;
            prevEngineeringStrain3D.mEngineeringStrain[2] -= deltaStrain;
*/
        }

		double deltaEqPlasticStrain;
		Eigen::Matrix<double,6,6> dSigmadEpsilon;
		Eigen::Matrix<double,6,6> dEpsilonPdEpsilon;
		Eigen::Matrix<double,6,1> newEpsilonP;
		Eigen::Matrix<double,6,1> newStress;

		bool uniaxialStressCondition(false);

		Eigen::Matrix<double,1,1> rdSigmadEpsilonRed; //reduced consistent tangent for the 1D problem
		while (uniaxialStressCondition==false)
		{
			// perform return mapping for the plasticity model
			NuTo::Error::eError error = this->ReturnMapping3D(
					engineeringStrain3D,
					prevEpsilonP,
					prevEngineeringStrain3D,
					newStress,
					newEpsilonP,
					deltaEqPlasticStrain,
					dSigmadEpsilon,
					dEpsilonPdEpsilon,
					rElement->GetStructure()->GetLogger());
			if (error!=Error::SUCCESSFUL)
				return error;

			//check, of uniaxial stress condition is reached
			double normStress(newStress.block<5,1>(1,0).norm());
			if (normStress<1e-10)
				uniaxialStressCondition=true;
			else
			{
				Eigen::Matrix<double,5,1>::Map(&(engineeringStrain3D.mEngineeringStrain[1]),5,1) +=
						dSigmadEpsilon.block<5,5>(1,1).inverse()*
						(newStress.block<5,1>(1,0)-
								dSigmadEpsilon.block<5,1>(1,0)*
								Eigen::Matrix<double,1,1>::Map(engineeringStrain3D.mEngineeringStrain,1,1));
			}
		}

		std::cout << "new plastic strain " << newEpsilonP << std::endl;

        //nonlocal damage
		if(rConstitutiveInput.find(NuTo::Constitutive::eInput::NONLOCAL_DAMAGE)==rConstitutiveInput.end())
            throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate1D] nonlocal damage needed to evaluate engineering strain1d.");
        double nonlocalDamage(rConstitutiveInput.find(NuTo::Constitutive::eInput::NONLOCAL_DAMAGE)->second->GetDamage());

        //new equivalent plastic strain
        double kappa = oldStaticData->mKappa+deltaEqPlasticStrain;
        //local damage
        double localDamage(1-exp(-kappa/mEpsilonF));

        //set this to true, if update is in the map, perform the update after all other outputs have been calculated
        bool performUpdateAtEnd(false);

        for (std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
                itOutput != rConstitutiveOutput.end(); itOutput++)
        {
            switch(itOutput->first)
            {
            case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_1D:
            {
                //calculate engineering stress
                EngineeringStress1D engineeringStress1D(itOutput->second->GetEngineeringStress1D());

                engineeringStress1D.mEngineeringStress = (1.-nonlocalDamage)*newStress(0);
            }
            break;
            case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_3D:
            {
                //this is for the visualize routines
                EngineeringStress3D& engineeringStress3D(itOutput->second->GetEngineeringStress3D());

                engineeringStress3D.mEngineeringStress[0] = (1.-nonlocalDamage)*newStress(0);
                engineeringStress3D.mEngineeringStress[1] = 0.;
                engineeringStress3D.mEngineeringStress[2] = 0.;
                engineeringStress3D.mEngineeringStress[3] = 0.;
                engineeringStress3D.mEngineeringStress[4] = 0.;
                engineeringStress3D.mEngineeringStress[5] = 0.;
            }
            break;
            case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D:
            {
        		Eigen::Matrix<double,1,1> dSigmadEpsilonRedInv
				 = dSigmadEpsilon.block<1,1>(0,0) +
				 dSigmadEpsilon.block<1,5>(0,1) * dSigmadEpsilon.block<5,5>(1,1).inverse() * dSigmadEpsilon.block<5,1>(1,0);
				ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
				tangent.SetSymmetry(true);
				double *localStiffData(tangent.mTangent);
				localStiffData[0] = (1.-nonlocalDamage) * dSigmadEpsilonRedInv(0,0);
            }
            break;
            case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_DAMAGE_1D:
            {
				ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
				tangent.SetSymmetry(true);
				double *localStiffData(tangent.mTangent);
				localStiffData[0] = -nonlocalDamage*newStress(0);
            }
            break;
            case NuTo::Constitutive::eOutput::D_NONLOCAL_DAMAGE_D_STRAIN_1D:
            {
				ConstitutiveTangentLocal<1,1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
				tangent.SetSymmetry(true);
				double *localStiffData(tangent.mTangent);
				localStiffData[0] = -nonlocalDamage*newStress(0);
            }
            break;
            case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_3D:
            {
                EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
    			engineeringStrain3D.mEngineeringStrain[0] = engineeringStrain3D.mEngineeringStrain[0];
    			engineeringStrain3D.mEngineeringStrain[1] = engineeringStrain3D.mEngineeringStrain[1];
    			engineeringStrain3D.mEngineeringStrain[2] = engineeringStrain3D.mEngineeringStrain[2];
    			engineeringStrain3D.mEngineeringStrain[3] = engineeringStrain3D.mEngineeringStrain[3];
    			engineeringStrain3D.mEngineeringStrain[4] = engineeringStrain3D.mEngineeringStrain[4];
    			engineeringStrain3D.mEngineeringStrain[5] = engineeringStrain3D.mEngineeringStrain[5];
            }
            break;
            case NuTo::Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_3D:
            {
            	EngineeringStrain3D& engineeringPlasticStrain(itOutput->second->GetEngineeringStrain3D());
            	engineeringPlasticStrain.mEngineeringStrain[0] = newEpsilonP(0);
            	engineeringPlasticStrain.mEngineeringStrain[1] = newEpsilonP(1);
            	engineeringPlasticStrain.mEngineeringStrain[2] = newEpsilonP(2);
            	engineeringPlasticStrain.mEngineeringStrain[3] = newEpsilonP(3);
            	engineeringPlasticStrain.mEngineeringStrain[4] = newEpsilonP(4);
            	engineeringPlasticStrain.mEngineeringStrain[5] = newEpsilonP(5);
            }
            break;
            case NuTo::Constitutive::eOutput::DAMAGE:
            {
                itOutput->second->GetDamage().SetDamage(localDamage);
            }
            break;
            case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
            {
           		throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate3D] tmp_static_data has to be updated without any other outputs, call it separately.");
            }
    		break;
            case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            {
            	performUpdateAtEnd = true;
            }
            break;
            default:
                throw MechanicsException(std::string("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate3D] output object)") +
                        NuTo::Constitutive::OutputToString(itOutput->first) +
                        std::string(" could not be calculated, check the allocated material law and the section behavior."));
            }
        }

        //update history variables
        if (performUpdateAtEnd)
        {
            double energy = oldStaticData->GetPrevTotalEnergy();
            //calculate delta total energy (sigma1+sigma2)/2*delta_strain
            const EngineeringStress1D& prevStress(oldStaticData->GetPrevStress());
            const EngineeringStrain1D prevStrain(oldStaticData->GetPrevStrain());
            energy+=0.5*(
                (newStress(0)+prevStress.mEngineeringStress)*(engineeringStrain3D.mEngineeringStrain[0]-prevStrain.mEngineeringStrain));
            oldStaticData->mPrevStrain = engineeringStrain1D;
            oldStaticData->mPrevSigma.mEngineeringStress = newStress(0);
            oldStaticData->SetPrevTotalEnergy(energy);

            //! @brief previous (after update) total strain component in radial direction (e22=e33)
            oldStaticData->mEpsilonTotRadial = engineeringStrain3D.mEngineeringStrain[1];

            //! @brief plastic strain (e11 and the radial component of the plastic strain e22=e33)
            oldStaticData->mEpsilonP[0] = newEpsilonP(0);
            oldStaticData->mEpsilonP[1] = newEpsilonP(1);

            // update the parts of the static data that are not related to the temporary updates (from the nonlocal calculation)
            oldStaticData->mKappa = kappa;
        }
        return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamagePlasticityEngineeringStress::Evaluate2D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
/*
	// get section information determining which input on the constitutive level should be used
    const SectionBase* section(rElement->GetSection());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }

    if (section->GetType()!=Section::PLANE_STRAIN)
    	throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate] only plane strain is implemented.");

    EngineeringStrain2D engineeringStrain;
    // calculate engineering strain
    if(rConstitutiveInput.find(NuTo::Constitutive::eInput::DEFORMATION_GRADIENT_2D)==rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate] deformation gradient 2d needed to evaluate engineering strain2d.");
    const DeformationGradient2D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::eInput::DEFORMATION_GRADIENT_2D)->second->GetDeformationGradient2D());
    deformationGradient.GetEngineeringStrain(engineeringStrain);

    //a recalculation of the plastic strain is not necessary, since this has been performed at the previous iteration with the update of the nonlocaltmpstatic data
    //Get previous ip_data
    ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain *oldStaticData = (rElement->GetStaticData(rIp))->AsNonlocalDamagePlasticity2DPlaneStrain();

    // subtract local plastic strain that has been calculated within the updatetmpstaticdata routine
    double elastStrain[4];
    elastStrain[0] = engineeringStrain.mEngineeringStrain[0] - oldStaticData->mTmpEpsilonP[0];
    elastStrain[1] = engineeringStrain.mEngineeringStrain[1] - oldStaticData->mTmpEpsilonP[1];
    elastStrain[2] = engineeringStrain.mEngineeringStrain[2] - oldStaticData->mTmpEpsilonP[2];
    elastStrain[3] =  - oldStaticData->mTmpEpsilonP[3];

    // subtract thermal strain
    if (section->GetInputConstitutiveIsTemperature())
    {
        std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::eInput::TEMPERATURE));
        if (itInput==rConstitutiveInput.end())
            throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate2D] temperature needed to evaluate thermal engineering strain2d.");
        double temperature(itInput->second->GetTemperature());
        double deltaStrain(mThermalExpansionCoefficient * temperature);
        EngineeringStrain2D elasticEngineeringStrain;
        elastStrain[0] -= deltaStrain;
        elastStrain[1] -= deltaStrain;
        elastStrain[3] -= deltaStrain;
    }

    // calculate coefficients of the material matrix
    double C11, C12, C33;
    this->CalculateCoefficients3D(C11, C12, C33);

    //calculate scaled equivalent plastic strain (scaled by the element length)
    bool unloading(true);
    double kappa = CalculateNonlocalEquivalentPlasticStrain(rElement, rIp, unloading);

    //determine kappaUnscaled, which is a scaling factor related to the fracture energy
    double kappaD = CalculateKappaD();

    //calculate engineering stress and damage
    EngineeringStress2D engineeringStress;
    double omega(0);
    double oneMinusOmega(1.);
    if (mDamage)
    {
        omega = CalculateDamage(kappa, kappaD);
        oneMinusOmega = 1.-omega;
        engineeringStress.mEngineeringrStress(0) = oneMinusOmega*(C11 * elastStrain[0] + C12*(elastStrain[1]+elastStrain[3]));
        engineeringStress.mEngineeringrStress(1) = oneMinusOmega*(C11 * elastStrain[1] + C12*(elastStrain[0]+elastStrain[3]));
        engineeringStress.mEngineeringrStress(2) = oneMinusOmega*(C33*elastStrain[2]);
    }
    else
    {
    	engineeringStress.mEngineeringrStress(0) = (C11 * elastStrain[0] + C12*(elastStrain[1]+elastStrain[3]));
    	engineeringStress.mEngineeringrStress(1) = (C11 * elastStrain[1] + C12*(elastStrain[0]+elastStrain[3]));
    	engineeringStress.mEngineeringrStress(2) = (C33*elastStrain[2]);
    }

    //set this to true, if update is in the map, perform the update after all other outputs have been calculated
    bool performUpdateAtEnd(false);

    for (std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
            itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch(itOutput->first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_2D:
        {
            // copy Engineering stress
            itOutput->second->GetEngineeringStress2D() = engineeringStress;

        }
        break;
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_3D:
        {
            //this is for the visualize routines
            EngineeringStress3D& engineeringStress3D(itOutput->second->GetEngineeringStress3D());

            engineeringStress3D.mEngineeringrStress(0) = engineeringStress.mEngineeringrStress(0);
            engineeringStress3D.mEngineeringrStress(1) = engineeringStress.mEngineeringrStress(1);
            engineeringStress3D.mEngineeringrStress(2) = oneMinusOmega*(C11 * elastStrain[3] + C12*(elastStrain[0]+elastStrain[1]));
            engineeringStress3D.mEngineeringrStress(3) = engineeringStress.mEngineeringrStress(2);
            engineeringStress3D.mEngineeringrStress(4) = 0.;
            engineeringStress3D.mEngineeringrStress(5) = 0.;
        }
        break;
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D:
        {
            if (mDamage)
            {
                //calculate damage parameter from the equivalente plastic strain
                //it is scaled related to the length (based on the direction of the first principal plastic strain)
                double dOmegadKappa;
                double oneMinusOmega = 1.-CalculateDerivativeDamage(kappa, kappaD, dOmegadKappa);
                assert(fabs(oneMinusOmega+omega-1.)<1e-10);
                //rLogger << "omega " << 1.-oneMinusOmega << "\n";
                if (unloading)
                {
                	ConstitutiveOutputBase* outputBase(itOutput->second);
                	outputBase->SetLocalSolution(true);
                	ConstitutiveTangentLocal<3,3>& tangent(outputBase->GetSubMatrix_3x3(0).AsConstitutiveTangentLocal_3x3());
                    //rLogger << "unloading local stiffness " << C11 << " "<< C12 << " " << C33 << " omega " << 1.-oneMinusOmega << "\n";
                    tangent.SetSymmetry(true);
                    double *localStiffData(tangent.mTangent);
                    localStiffData[0] = oneMinusOmega * C11;
                    localStiffData[1] = oneMinusOmega * C12;
                    localStiffData[2] = 0.;

                    localStiffData[3] = localStiffData[1];
                    localStiffData[4] = localStiffData[0];
                    localStiffData[5] = 0.;

                    localStiffData[6] = 0.;
                    localStiffData[7] = 0.;
                    localStiffData[8] = oneMinusOmega * C33;
                }
                else
                {
                    //rLogger << "nonlocal stiffness, omega " << 1.-oneMinusOmega << "\n";
                	ConstitutiveOutputBase* outputBase(itOutput->second);
                	outputBase->SetLocalSolution(false);

                	//get nonlocal elements
                    const std::vector<const NuTo::ElementBase*>& nonlocalElements(rElement->GetNonlocalElements());

                    //calculate nonlocal plastic strain for the direction of the equivalent length
                    double nonlocalPlasticStrain[4];
                    CalculateNonlocalPlasticStrain(rElement, rIp, nonlocalPlasticStrain);

                    // calculate Engineering stress
                    Eigen::Matrix<double,3,1> sigmaElast;
                    sigmaElast(0) = C11 * elastStrain[0] + C12*(elastStrain[1]+elastStrain[3]);
                    sigmaElast(1) = C11 * elastStrain[1] + C12*(elastStrain[0]+elastStrain[3]);
                    sigmaElast(2) = C33 * elastStrain[2] ;

                    Eigen::Matrix<double,4,1> dKappadEpsilonP;
                    //calculate derivative of damage parameter with respect to local strains of all nonlocal integration points
                    int totalNonlocalIp(0);
                    for (int countNonlocalElement=0; countNonlocalElement<(int)nonlocalElements.size(); countNonlocalElement++)
                    {
                        const std::vector<double>& weights(rElement->GetNonlocalWeights(rIp,countNonlocalElement));

                        //rLogger << weights.size() << " "<< nonlocalElements[countNonlocalElement]->GetNumIntegrationPoints() << "\n";
                        assert((int)weights.size()==nonlocalElements[countNonlocalElement]->GetNumIntegrationPoints());

                        //Go through all the integration points
                        for (int theNonlocalIP=0; theNonlocalIP<(int)weights.size(); theNonlocalIP++, totalNonlocalIp++)
                        {
                            if (weights[theNonlocalIP]==0.)
                                continue;
                            assert(totalNonlocalIp<outputBase->GetNumSubMatrices());
                            ConstitutiveTangentLocal<3,3>& tangent(outputBase->GetSubMatrix_3x3(totalNonlocalIp));
                            double *localStiffData = tangent.mTangent;

                            //here the nonlocal delta is relevant
                            const ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain *NonlocalOldStaticData = (nonlocalElements[countNonlocalElement]->GetStaticData(theNonlocalIP))->AsNonlocalDamagePlasticity2DPlaneStrain();

                            double deltaEpsilonPxx = NonlocalOldStaticData->mTmpEpsilonP[0] - NonlocalOldStaticData->mEpsilonP[0];
                            double deltaEpsilonPyy = NonlocalOldStaticData->mTmpEpsilonP[1] - NonlocalOldStaticData->mEpsilonP[1];
                            double deltaEpsilonPxy = NonlocalOldStaticData->mTmpEpsilonP[2] - NonlocalOldStaticData->mEpsilonP[2];
                            double deltaEpsilonPzz = NonlocalOldStaticData->mTmpEpsilonP[3] - NonlocalOldStaticData->mEpsilonP[3];

                            double deltaEpsilonPEq = sqrt(deltaEpsilonPxx*deltaEpsilonPxx+deltaEpsilonPyy*deltaEpsilonPyy+
                                                            0.5*deltaEpsilonPxy*deltaEpsilonPxy+deltaEpsilonPzz*deltaEpsilonPzz);
                            if (deltaEpsilonPEq>0)
                            {
                                double factor(NonlocalOldStaticData->mTmpLeq/deltaEpsilonPEq);
                                dKappadEpsilonP(0) = factor*deltaEpsilonPxx;
                                dKappadEpsilonP(1) = factor*deltaEpsilonPyy;
                                dKappadEpsilonP(2) = factor*0.5*deltaEpsilonPxy;
                                dKappadEpsilonP(3) = factor*deltaEpsilonPzz;
                            }
                            else
                            {
                                double factor(NonlocalOldStaticData->mTmpLeq);
                                dKappadEpsilonP(0) = factor;
                                dKappadEpsilonP(1) = factor;
                                dKappadEpsilonP(2) = factor*0.5;
                                dKappadEpsilonP(3) = factor;
                            }

                            //calculate dOmegadepsilon, instead of using transpose, just declare a RowMajor storage
                            Eigen::Matrix<double,3,1> minusdOmegadEpsilon (Eigen::Matrix<double,4,4,Eigen::RowMajor>::Map(NonlocalOldStaticData->mTmpdEpsilonPdEpsilon,4,4).topLeftCorner<3,4>() * dKappadEpsilonP);

                            //the second part includes the dependence of leq as a function of the local stress at the nonlocal integration point
                            minusdOmegadEpsilon += Eigen::Matrix<double,3,1>::Map(NonlocalOldStaticData->mTmpdLeqdEpsilon,3)*deltaEpsilonPEq;

                            minusdOmegadEpsilon *= -dOmegadKappa * weights[theNonlocalIP];
                            //rLogger<< "dOmegadepsilon/weight analytic" << "\n" << -minusdOmegadEpsilon/weights[theNonlocalIP] << "\n";
                            //rLogger << "weight " << weights[theNonlocalIP] << "\n";

                            Eigen::Matrix<double,3,Eigen::Dynamic>::Map(localStiffData,3, 3) = sigmaElast * minusdOmegadEpsilon.transpose() ;

                            if (nonlocalElements[countNonlocalElement]==rElement && rIp == theNonlocalIP)
                            {

                                const double *tmpdEP(oldStaticData->mTmpdEpsilonPdEpsilon);
                                localStiffData[0] += oneMinusOmega*(C11*(1.-tmpdEP[0]) - C12*(tmpdEP[1] + tmpdEP[3]));
                                localStiffData[1] += oneMinusOmega*(C12*(1.-tmpdEP[0] - tmpdEP[3]) - C11*tmpdEP[1]);
                                localStiffData[2] += oneMinusOmega*(-C33*tmpdEP[2]);
                                localStiffData[3] += oneMinusOmega*(-C11*tmpdEP[4] + C12*(1.-tmpdEP[5]- tmpdEP[7]));
                                localStiffData[4] += oneMinusOmega*(-C12*(tmpdEP[4] + tmpdEP[7]) + C11*(1.-tmpdEP[5]));
                                localStiffData[5] += oneMinusOmega*(-C33*tmpdEP[6]);
                                localStiffData[6] += oneMinusOmega*(-C11*tmpdEP[8] - C12*(tmpdEP[9]+tmpdEP[11]));
                                localStiffData[7] += oneMinusOmega*(-C12*(tmpdEP[8]+ tmpdEP[11]) - C11*tmpdEP[9]);
                                localStiffData[8] += oneMinusOmega*(C33*(1.-tmpdEP[10]));
                            }
                        }
                    }
                }
            }
            else
            {
            	ConstitutiveOutputBase* outputBase(itOutput->second);
            	outputBase->SetLocalSolution(true);
            	ConstitutiveTangentLocal<3,3>& tangent(outputBase->GetSubMatrix_3x3(0).AsConstitutiveTangentLocal_3x3());
                tangent.SetSymmetry(true);
                double *localStiffData(tangent.mTangent);

                const double *tmpdEP(oldStaticData->mTmpdEpsilonPdEpsilon);
                localStiffData[0] = C11*(1.-tmpdEP[0]) - C12*(tmpdEP[1] + tmpdEP[3]);
                localStiffData[1] = C12*(1.-tmpdEP[0] - tmpdEP[3]) - C11*tmpdEP[1];
                localStiffData[2] = -C33*tmpdEP[2];
                localStiffData[3] = -C11*tmpdEP[4] + C12*(1.-tmpdEP[5]- tmpdEP[7]);
                localStiffData[4] = -C12*(tmpdEP[4] + tmpdEP[7]) + C11*(1.-tmpdEP[5]);
                localStiffData[5] = -C33*tmpdEP[6];
                localStiffData[6] = -C11*tmpdEP[8] - C12*(tmpdEP[9]+tmpdEP[11]);
                localStiffData[7] = -C12*(tmpdEP[8]+ tmpdEP[11]) - C11*tmpdEP[9];
                localStiffData[8] = C33*(1.-tmpdEP[10]);
            }
        }
        break;
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_3D:
        {
            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
			engineeringStrain3D.mEngineeringStrain[0] = engineeringStrain.mEngineeringStrain[0];
			engineeringStrain3D.mEngineeringStrain[1] = engineeringStrain.mEngineeringStrain[1];
			engineeringStrain3D.mEngineeringStrain[2] = 0;
			engineeringStrain3D.mEngineeringStrain[3] = engineeringStrain.mEngineeringStrain[2];
			engineeringStrain3D.mEngineeringStrain[4] = 0.;
			engineeringStrain3D.mEngineeringStrain[5] = 0.;
        }
        break;
        case NuTo::Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_3D:
        {
        	EngineeringStrain3D& engineeringPlasticStrain(itOutput->second->GetEngineeringStrain3D());
        	engineeringPlasticStrain.mEngineeringStrain[0] = oldStaticData->mTmpEpsilonP[0];
        	engineeringPlasticStrain.mEngineeringStrain[1] = oldStaticData->mTmpEpsilonP[1];
        	engineeringPlasticStrain.mEngineeringStrain[2] = oldStaticData->mTmpEpsilonP[3];
        	engineeringPlasticStrain.mEngineeringStrain[3] = oldStaticData->mTmpEpsilonP[2];
        	engineeringPlasticStrain.mEngineeringStrain[4] = 0.;
        	engineeringPlasticStrain.mEngineeringStrain[5] = 0.;
        }
        break;
        case NuTo::Constitutive::eOutput::DAMAGE:
        {
            itOutput->second->GetDamage().SetDamage(omega);        }
        break;
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        {
        	if (rConstitutiveOutput.size()!=1)
        		throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate3D] tmp_static_data has to be updated without any other outputs, call it separately.");

			double rDeltaEqPlasticStrain;
			Eigen::Matrix<double,4,4> rdEpsilonPdEpsilon;
			Eigen::Matrix<double,4,1> rNewEpsilonP;
			Eigen::Matrix<double,4,1> rNewStress;

			// perform return mapping for the plasticity model
			NuTo::Error::eError error = this->ReturnMapping2D(
					engineeringStrain,
					oldStaticData->mEpsilonP,
					oldStaticData->mPrevStrain,
					rNewStress,
					rNewEpsilonP,
					rDeltaEqPlasticStrain,
					rdEpsilonPdEpsilon,
					rElement->GetStructure()->GetLogger());

			if (error!=Error::SUCCESSFUL)
				return error;

			std::cout << "new plastic strain " << rNewEpsilonP << std::endl;

			// update the temporary parts of the static data
			Eigen::Matrix<double,4,1>::Map(oldStaticData->mTmpEpsilonP,4) = rNewEpsilonP;

			//determine equivalente length in zz and plane direction
			//using the local elastic stress
			Eigen::Matrix<double,4,1> dLdSigma;
			oldStaticData->mTmpLeq = CalculateDerivativeEquivalentLength2D(rElement,rNewStress,dLdSigma);

			oldStaticData->mTmpKappa = oldStaticData->mKappa + rDeltaEqPlasticStrain*oldStaticData->mTmpLeq;

			//rLogger << "tmpKappa " << oldStaticData->mTmpKappa << "\n";
			Eigen::Matrix<double,4,4>::Map(oldStaticData->mTmpdEpsilonPdEpsilon,4,4) = rdEpsilonPdEpsilon;

			// calculate coefficients of the linear elastic material matrix
			Eigen::Matrix<double,4,4> ElasticStiffness;
			ElasticStiffness << C11, C12, 0.,  C12,
								C12, C11, 0.,  C12,
								 0.,  0., C33, 0.,
								C12, C12, 0.,  C11;

			// calculate derivative of eq length with respect to local strain
			Eigen::Matrix<double,4,1>::Map(oldStaticData->mTmpdLeqdEpsilon,4) = dLdSigma.transpose() * ElasticStiffness * (Eigen::Matrix<double,4,4>::Identity() - Eigen::Matrix<double,4,4>::Map(oldStaticData->mTmpdEpsilonPdEpsilon,4,4));
        }
		break;
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
        	performUpdateAtEnd = true;
        }
        break;
        default:
            throw MechanicsException(std::string("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate3D] output object)") +
                    NuTo::Constitutive::OutputToString(itOutput->first) +
                    std::string(" culd not be calculated, check the allocated material law and the section behavior."));
        }
    }

    //update history variables
    if (performUpdateAtEnd)
    {
        double energy = oldStaticData->GetPrevTotalEnergy();
        //calculate delta total energy (sigma1+sigma2)/2*delta_strain
        const EngineeringStress2D& prevStress(oldStaticData->GetPrevStress());
        const EngineeringStrain2D prevStrain(oldStaticData->GetPrevStrain());
        energy+=0.5*(
            (engineeringStress.mEngineeringrStress(0)+prevStress.mEngineeringrStress(0))*(engineeringStrain.mEngineeringStrain[0]-prevStrain.mEngineeringStrain[0])+
            (engineeringStress.mEngineeringrStress(1)+prevStress.mEngineeringrStress(1))*(engineeringStrain.mEngineeringStrain[1]-prevStrain.mEngineeringStrain[1])+
            (engineeringStress.mEngineeringrStress(2)+prevStress.mEngineeringrStress(2))*(engineeringStrain.mEngineeringStrain[2]-prevStrain.mEngineeringStrain[2])
            );
        oldStaticData->mPrevStrain = engineeringStrain;
        oldStaticData->mPrevSigma = engineeringStress;
        oldStaticData->SetPrevTotalEnergy(energy);

        // update the parts of the static data that are not related to the temporary updates (from the nonlocal calculation)
        oldStaticData->mKappa = oldStaticData->mTmpKappa;

//        Eigen::Matrix<double,4,1>::Map(oldStaticData->mEpsilonP,4,1) = Eigen::Matrix<double,4,1>::Map(oldStaticData->mTmpEpsilonP,4,1);;
    }
*/
	throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate2D] not implemented for 2D.");
    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamagePlasticityEngineeringStress::Evaluate3D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::Evaluate3D] not implemented for 3D.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
        const ElementBase* rElement) const
{
	return new ConstitutiveStaticDataGradientDamagePlasticity1D();
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
        const ElementBase* rElement) const
{
/*    if (rElement->GetSection()==0)
        throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Section required to distinguish between plane stress and plane strain and thickness information.");
    if (rElement->GetSection()->GetType()==NuTo::Section::PLANE_STRESS)
        throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for plane stress.");
    else
        return new ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain();
*/
    throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D] Nonlocal damage plasticity model not implemented for 2D.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
        const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D] Nonlocal damage plasticity model not implemented for 3D.");
    //return new ConstitutiveStaticDataGradientDamagePlasticity3D();
}


// calculate coefficients of the material matrix
void NuTo::GradientDamagePlasticityEngineeringStress::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE/((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE/(2*(1.0 + this->mNu));
}

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get density
//! @return ... density
double NuTo::GradientDamagePlasticityEngineeringStress::GetDensity() const
{
    return this->mRho;
}

//! @brief ... set density
//! @param rRho ... density
void NuTo::GradientDamagePlasticityEngineeringStress::SetDensity(double rRho)
{
    this->CheckDensity(rRho);
    this->mRho = rRho;
    this->SetParametersValid();
}

//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::GradientDamagePlasticityEngineeringStress::GetYoungsModulus() const
{
    return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::GradientDamagePlasticityEngineeringStress::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::GradientDamagePlasticityEngineeringStress::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::GradientDamagePlasticityEngineeringStress::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}

//! @brief ... get nonlocal radius
//! @return ... nonlocal radius
double NuTo::GradientDamagePlasticityEngineeringStress::GetNonlocalRadius() const
{
    return mNonlocalRadius;
}

//! @brief ... set nonlocal radius
//! @param rRadius...  nonlocal radius
void NuTo::GradientDamagePlasticityEngineeringStress::SetNonlocalRadius(double rNonlocalRadius)
{
    this->CheckNonlocalRadius(rNonlocalRadius);
    this->mNonlocalRadius = rNonlocalRadius;
    this->SetParametersValid();
}
//! @brief ... get tensile strength
//! @return ... tensile strength
double NuTo::GradientDamagePlasticityEngineeringStress::GetTensileStrength() const
{
    return mTensileStrength;
}

//! @brief ... set tensile strength
//! @param rTensileStrength...  tensile strength
void NuTo::GradientDamagePlasticityEngineeringStress::SetTensileStrength(double rTensileStrength)
{
    this->CheckNonlocalRadius(rTensileStrength);
    this->mTensileStrength = rTensileStrength;
    this->SetParametersValid();
}

//! @brief ... get compressive strength
//! @return ... compressive strength
double NuTo::GradientDamagePlasticityEngineeringStress::GetCompressiveStrength() const
{
    return mCompressiveStrength;
}

//! @brief ... set compressive strength
//! @param rCompressiveStrength...  compressive strength
void NuTo::GradientDamagePlasticityEngineeringStress::SetCompressiveStrength(double rCompressiveStrength)
{
    this->CheckNonlocalRadius(rCompressiveStrength);
    this->mCompressiveStrength = rCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get biaxial compressive strength
//! @return ... biaxial compressive strength
double NuTo::GradientDamagePlasticityEngineeringStress::GetBiaxialCompressiveStrength() const
{
    return mBiaxialCompressiveStrength;
}

//! @brief ... set biaxial compressive strength
//! @param rBiaxialCompressiveStrength...  biaxial compressive strength
void NuTo::GradientDamagePlasticityEngineeringStress::SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength)
{
    this->CheckNonlocalRadius(rBiaxialCompressiveStrength);
    this->mBiaxialCompressiveStrength = rBiaxialCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get fracture energy
//! @return ... fracture energy
double NuTo::GradientDamagePlasticityEngineeringStress::GetFractureEnergy() const
{
    return mFractureEnergy;
}

//! @brief ... set fracture energy
//! @param rFractureEnergy... fracture energy
void NuTo::GradientDamagePlasticityEngineeringStress::SetFractureEnergy(double rFractureEnergy)
{
    this->CheckFractureEnergy(rFractureEnergy);
    this->mFractureEnergy = rFractureEnergy;
    this->SetParametersValid();
}

//! @brief ... get thermal expansion coefficient
//! @return ... thermal expansion coefficient
double NuTo::GradientDamagePlasticityEngineeringStress::GetThermalExpansionCoefficient() const
{
    return mThermalExpansionCoefficient;
}

//! @brief ... set thermal expansion coefficient
//! @param rNu ... thermal expansion coefficient
void NuTo::GradientDamagePlasticityEngineeringStress::SetThermalExpansionCoefficient(double rAlpha)
{
    this->CheckThermalExpansionCoefficient(rAlpha);
    this->mThermalExpansionCoefficient = rAlpha;
    this->SetParametersValid();
}
///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::GradientDamagePlasticityEngineeringStress::GetType() const
{
    return NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::GradientDamagePlasticityEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return false;
    case NuTo::Element::PLANE2D3N:
        return false;
    case NuTo::Element::PLANE2D4N:
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
        return true;
    default:
        return false;
    }
}

//! @brief ... check if density is non negative
//! @param rE ... Young's modulus
void NuTo::GradientDamagePlasticityEngineeringStress::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckDensity] The density must not be negative.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::GradientDamagePlasticityEngineeringStress::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::GradientDamagePlasticityEngineeringStress::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check if the nonlocal radius is positive
//! @param rRadius ... nonlocal radius
void NuTo::GradientDamagePlasticityEngineeringStress::CheckNonlocalRadius(double rRadius) const
{
    if (rRadius <= 0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckNonlocalRadius] Nonlocal radius must be positive.");
    }
}
//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... nonlocal radius
void NuTo::GradientDamagePlasticityEngineeringStress::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if compressive strength is positive
//! @param rRadius ... compressive strength
void NuTo::GradientDamagePlasticityEngineeringStress::CheckCompressiveStrength(double rCompressiveStrength) const
{
    if (rCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckCompressiveStrength] The compressive strength must be a positive value.");
    }
}

//! @brief ... check if biaxial compressive strength is positive
//! @param rBiaxialCompressiveStrength ... biaxial compressive strength
void NuTo::GradientDamagePlasticityEngineeringStress::CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const
{
    if (rBiaxialCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be a positive value.");
    }
    if (rBiaxialCompressiveStrength <= mCompressiveStrength)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be higher than the uniaxial compressive strength.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::GradientDamagePlasticityEngineeringStress::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}

//! @brief ... check thermal expansion coefficient
//! @param rAlpha ... thermal expansion coefficient
void NuTo::GradientDamagePlasticityEngineeringStress::CheckThermalExpansionCoefficient(double rAlpha) const
{
}

#define ACTIVE true
#define INACTIVE false
#define toleranceResidual 1e-7      //tolerance to decide whether the Newton iteration has converged
#define toleranceYieldSurface 1e-9  //tolerance whether a point is on the yield surface or not (multiplied by the tensile strength)
#define toleranceDeterminant 1e-50  //tolerance to decide if a matrix is not invertible (only required in the debugging version, be careful here with the units)
#define maxSteps 25                 //maximum number of Newton iterations, until it is decided that there is no convergence and a cutback is performed
#define minCutbackFactor 1e-3       //minimum cutback factor for the application of the total strain in steps
#define minCutbackFactorLS 2e-3     //minimum cutback factor used for the linesearch in the Newton iteration


//! @brief ... performs the return mapping procedure for the plasticity model
//! @param rStrain              ... current total strain
//! @param rPrevPlasticStrain   ... previous plastic strain (history variable)
//! @param rPrevTotalStrain     ... previous total strain (history variable)
//! @param rPrevEqPlasticStrain ... previous equiavalente plastic strain (history variable)
//! @param rEpsilonP            ... new plastic strain after return mapping
//! @param rEqPlasticStrain     ... new equivalente olastic strain after return mapping
//! @param rdEpsilonPdEpsilon   ... new derivative of current plastic strain with respect to the total strain
NuTo::Error::eError NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping1D(
        const EngineeringStrain1D& rStrain,
        const double rPrevPlasticStrain[2],
        const EngineeringStrain1D& rPrevTotalStrain,
        Eigen::Matrix<double,1,1>& rStress,
        Eigen::Matrix<double,2,1>& rEpsilonP,
        double& rDeltaEqPlasticStrain,
        Eigen::Matrix<double,1,1>& rdEpsilonPdEpsilon,
        NuTo::Logger& rLogger)const
{
    return Error::SUCCESSFUL;
}

//! @brief ... performs the return mapping procedure for the plasticity model
//! @param rStrain              ... current total strain
//! @param rPrevPlasticStrain   ... previous plastic strain (history variable)
//! @param rPrevTotalStrain     ... previous total strain (history variable)
//! @param rPrevEqPlasticStrain ... previous equiavalente plastic strain (history variable)
//! @param rEpsilonP            ... new plastic strain after return mapping
//! @param rEqPlasticStrain     ... new equivalente olastic strain after return mapping
//! @param rdEpsilonPdEpsilon   ... new derivative of current plastic strain with respect to the total strain
NuTo::Error::eError NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D(
        const EngineeringStrain2D& rStrain,
        const double rPrevPlasticStrain[4],
        const EngineeringStrain2D& rPrevTotalStrain,
        Eigen::Matrix<double,4,1>& rStress,
        Eigen::Matrix<double,4,1>& rEpsilonP,
        double& rDeltaEqPlasticStrain,
        Eigen::Matrix<double,4,4>& rdEpsilonPdEpsilon,
        NuTo::Logger& rLogger)const
{
/*
	double e_mod = mE; //modify that one in the case of random fields
    double   nu  = mNu;
    double f_ct  = mTensileStrength;
    double f_c1  = mCompressiveStrength;
    double f_c2  = mBiaxialCompressiveStrength;

    assert(f_c2>f_c1);
    assert(f_c1>0);
    assert(f_c2>0);
    assert(e_mod>0);

    // ******************************************************************
    // *    F_BETA:    required by DRUCKER-PRAGER yield surface               *
    // ******************************************************************
    double BETA = sqrt3*(f_c2-f_c1) / (2*f_c2-f_c1);
    double H_P  = f_c2*f_c1 / (sqrt3*(2*f_c2-f_c1));

    //! @brief strain currently solved for the plastic strains, in general equal to  rStrain, but if applied in steps, it's smaller
    Eigen::Vector4d curTotalStrain;
    //! @brief previous plastic strain, either from previous equilibrium (static data) or if applied in steps, previous converged state
    Eigen::Vector4d lastPlastStrain;
    //! @brief plastic strain in the line search
    Eigen::Vector4d epsilonPLS;
    //! @brief previous eq plastic strain, either from previous equilibrium (static data) or if applied in steps, previous converged state
    double lastDeltaEqPlasticStrain;
    //! @brief residual in the return mapping procedure
    Eigen::Vector4d residual;
    //! @brief residual in the return mapping procedure within linesearch
    Eigen::Vector4d residualLS;
    //! @brief full stress increment within one iteration of the return mapping, might be applied in steps in the succeeding line search
    Eigen::Vector4d deltaStress;
    //! @brief total strain increment between strain from previous static data and new total strain
    Eigen::Vector4d deltaStrain;
    //! @brief plastic strain increment within the linesearch
    Eigen::Vector4d deltaPlasticStrainLS;
    //! @brief plastic strain increment between to load steps (internal load splitting)
    Eigen::Vector4d deltaPlasticStrain;
    //! @brief trial stress of the first iteration
    Eigen::Vector4d initTrialStress;
    //! @brief trial stress in the line search
    Eigen::Vector4d stressLS;
    //! @brief trial stress
    Eigen::Vector4d trialStress;
    //! @brief elastic strain
    Eigen::Vector4d elasticStrain;
    //! @brief elastic strain in line search
    Eigen::Vector4d elasticStrainLS;
    //! @brief plastic multiplier
    Eigen::Matrix<double,2,1> deltaGamma;
    //! @brief plastic multiplier
    Eigen::Matrix<double,2,1> deltaGammaLS;
    //! @brief increment of plastic multiplier in return mapping procedure
    Eigen::Matrix<double,2,1> delta2Gamma;
    //! @brief yield condition
    Eigen::Matrix<double,2,1> yieldCondition;
    //! @brief yield condition in line search
    Eigen::Matrix<double,2,1> yieldConditionLS(0,0);
    //! @brief yield condition at the first iteration
    Eigen::Matrix<double,2,1> initYieldCondition;
    //! @brief flag that indicates if a yield function is active or not
    Eigen::Matrix<bool,2,1> yieldConditionFlag;
    //! @brief (dF/dsigma)T * Hessian * dF/dsigma
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matG;
    //! @brief ((dF/dsigma)T * Hessian * dF/dsigma )^-1
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matGInv;
    //! @brief algorithmic modulus DElastInv+delta_gamm*d2F/d2Sigma
    Eigen::Matrix<double,4,4> hessian;
    //! @brief temporary matrix
    Eigen::Matrix<double,4,4> tmpStiffness;
    //! @brief first derivatives of the yield functions
    std::vector<Eigen::Matrix<double,4,1> > dF_dsigma;
    //! @brief second derivatives of the yield functions
    std::vector<Eigen::Matrix<double,4,4> > d2F_d2Sigma;
    //! @brief algorithmic modulus * dF_dsigma
    std::vector<Eigen::Matrix<double,4,1> > vectorN;
    //! @brief number of active yield functions
    int numActiveYieldFunctions;

    // for the application of strains in steps, calculate the total strain increment to be applied
    deltaStrain(0) = rStrain.mEngineeringStrain[0]-rPrevTotalStrain.mEngineeringStrain[0];
    deltaStrain(1) = rStrain.mEngineeringStrain[1]-rPrevTotalStrain.mEngineeringStrain[1];
    deltaStrain(2) = rStrain.mEngineeringStrain[2]-rPrevTotalStrain.mEngineeringStrain[2];
    deltaStrain(3) = 0;

    // initialize last plastic strain and last converged stress
    lastPlastStrain << rPrevPlasticStrain[0] , rPrevPlasticStrain[1] ,rPrevPlasticStrain[2] ,rPrevPlasticStrain[3];
    rDeltaEqPlasticStrain = 0;
    // *****************************************************************
    //                   elastic matrix generation                     *
    // *****************************************************************
    //! @brief elastic stiffness
    Eigen::Matrix<double,4,4> dElast;
    //! @brief inverse elastic stiffness
    Eigen::Matrix<double,4,4> dElastInv;

    {
        double factor = e_mod/((1.+nu)*(1.-2.*nu));
        double oneminusnufactor = (1-nu)*factor;
        double nufactor = nu*factor;

        dElast <<  oneminusnufactor , nufactor         , 0.              , nufactor ,
                   nufactor         , oneminusnufactor , 0.              , nufactor ,
                   0.               , 0.               , (0.5-nu)*factor , 0.,
                   nufactor         , nufactor         , 0.              ,oneminusnufactor;

        factor = 1./e_mod;
        double minusnufactor = -nu*factor;
        dElastInv   << factor        , minusnufactor , 0.               , minusnufactor,
                       minusnufactor , factor        , 0.               , minusnufactor,
                       0.            , 0.            , 2.*factor*(1.+nu), 0.,
                       minusnufactor , minusnufactor , 0                , factor;

    }
    //! @brief delta load factor for the previous iteration
    double deltaCutbackFactorExternal(1.);

    //! @brief current load factor (between 0 and 1) to apply the total strain increment in steps
    double cutbackFactorExternal(deltaCutbackFactorExternal);

    //! @brief flag to determine if the iteration is finished (converged at  cutbackFactorExternal=1)
    bool convergedExternal(false);

    int numberOfExternalCutbacks(0);
    int numberOfInternalIterations(0);
#ifdef ENABLE_DEBUG
    int prevNumberOfInternalIterations(0);
#endif
    lastDeltaEqPlasticStrain = 0.;
    while (cutbackFactorExternal>minCutbackFactor && !convergedExternal)
    {
        numberOfExternalCutbacks++;

        curTotalStrain(0) = rPrevTotalStrain.mEngineeringStrain[0]+cutbackFactorExternal*deltaStrain(0);
        curTotalStrain(1) = rPrevTotalStrain.mEngineeringStrain[1]+cutbackFactorExternal*deltaStrain(1);
        curTotalStrain(2) = rPrevTotalStrain.mEngineeringStrain[2]+cutbackFactorExternal*deltaStrain(2);
        curTotalStrain(3) = 0.;

#ifdef ENABLE_DEBUG
        rLogger << "\n" << "curTotalStrain" << curTotalStrain.transpose() << "\n" << "\n";
#endif

        // checks the convergence of the Newton iteration for a prescribed current strain
        bool convergedInternal(false);
        try
        {
            //resize yield condition vector
            int numYieldSurfaces=2;
            yieldCondition.setZero(numYieldSurfaces);

            //elastic strain, stress and d_matrix
            elasticStrain = curTotalStrain - lastPlastStrain;

            //TODO just use the upper part for new EigenVersion 3.0
            //trialStress = dElast.selfadjointView<Eigen::Upper>()*elasticStrain;
            initTrialStress = dElast*elasticStrain;
#ifdef ENABLE_DEBUG
            rLogger << "initTrialStress " << "\n" << initTrialStress.transpose() << "\n" << "\n";
#endif

            //calculate yield condition
            //Drucker Prager
            initYieldCondition(0) = YieldSurfaceDruckerPrager2D(initTrialStress, BETA, H_P);

            //rounded Rankine
            initYieldCondition(1) = YieldSurfaceRankine2DRounded(initTrialStress, f_ct);

#ifdef ENABLE_DEBUG
            rLogger << "initYieldCondition " << "\n" << initYieldCondition(0) << " " << initYieldCondition(1) << "\n";
#endif


            if (initYieldCondition(0)<-toleranceYieldSurface*f_ct && initYieldCondition(1)<-toleranceYieldSurface*f_ct)
            {
                // *************************************************
                // *  thus we have elastic -------------> elastic  *
                // *************************************************
#ifdef ENABLE_DEBUG
                rLogger << "linear elastic step" << "\n" << "\n";
#endif
                convergedInternal = true;
                rEpsilonP =  lastPlastStrain;
                rDeltaEqPlasticStrain = 0;
                trialStress = initTrialStress;
                if (cutbackFactorExternal==1)
                {
                    rdEpsilonPdEpsilon.setZero(4,4);
                    rStress = trialStress;
                    return Error::SUCCESSFUL;
;
                }
            }
            else
            {
#ifdef ENABLE_DEBUG
                rLogger << "plastic step" << "\n" << "\n";
#endif
            }

            // perform return mapping
            dF_dsigma.resize(numYieldSurfaces);
            d2F_d2Sigma.resize(numYieldSurfaces);

            // initialize plastic multiplier
            deltaGamma.setZero(numYieldSurfaces);
            delta2Gamma.setZero(numYieldSurfaces);
            yieldConditionFlag.setZero(numYieldSurfaces);

            for (int fixedYieldConditions=0; fixedYieldConditions<3 && convergedInternal==false; fixedYieldConditions++)
            {
                switch(fixedYieldConditions)
                {
                case 0:
#ifdef ENABLE_DEBUG
                    rLogger<< "pure Rankine" << "\n";
#endif
                    deltaGamma(1) = 0;
                    if (initYieldCondition(1)<-toleranceYieldSurface*f_ct)
                        continue;
                    yieldConditionFlag(0) = INACTIVE;
                    yieldConditionFlag(1) = ACTIVE;
                    numActiveYieldFunctions = 1;
                    break;
                case 1:
#ifdef ENABLE_DEBUG
                    rLogger<< "combined" << "\n";
#endif
                    if (initYieldCondition(0)<-toleranceYieldSurface*f_ct || initYieldCondition(1)<-toleranceYieldSurface*f_ct)
                        continue;
                    yieldConditionFlag(0) = ACTIVE;
                    yieldConditionFlag(1) = ACTIVE;
                    deltaGamma(0) = 0;
                    deltaGamma(1) = 0;
                    numActiveYieldFunctions = 2;
                    break;
                case 2:
#ifdef ENABLE_DEBUG
                    rLogger<< "pure DP" << "\n";
#endif
                    if (initYieldCondition(0)<-toleranceYieldSurface*f_ct)
                        continue;
                    yieldConditionFlag(0) = ACTIVE;
                    yieldConditionFlag(1) = INACTIVE;
                    deltaGamma(0) = 0;
                    numActiveYieldFunctions = 1;
                    break;
                default:
                    throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] programming error - should not happen.");
                }

                for (int iteration = 0; iteration < maxSteps && convergedInternal==false; iteration++)
                {
                    numberOfInternalIterations++;
                    if (iteration==0)
                    {
                        rEpsilonP = lastPlastStrain;
                        rDeltaEqPlasticStrain = lastDeltaEqPlasticStrain;
                        yieldCondition= initYieldCondition;
                        trialStress = initTrialStress;
                    }
                    else
                    {
                        //trial stress is the last stress state from the previous line search
                        if (yieldConditionFlag(0)==ACTIVE)
                            yieldCondition(0) = yieldConditionLS(0);
                        else
                            yieldCondition(0) = YieldSurfaceDruckerPrager2D(trialStress, BETA, H_P);

                        if (yieldConditionFlag(1)==ACTIVE)
                            yieldCondition(1) = yieldConditionLS(1);
                        else
                            yieldCondition(1) = YieldSurfaceRankine2DRounded(trialStress, f_ct);

                    }
#ifdef ENABLE_DEBUG
                    rLogger << "trialStress " <<  "\n" << trialStress.transpose() << "\n" << "\n";
                    rLogger << "yieldCondition " <<  "\n" << yieldCondition.transpose() << "\n" << "\n";
#endif

                    // DP
                    if (yieldConditionFlag(0)==ACTIVE)
                    {
                        if (!YieldSurfaceDruckerPrager2DDerivatives(dF_dsigma[0],&(d2F_d2Sigma[0]),trialStress,BETA))
                        {
                            //no convergence, decrease line search step
                            iteration = maxSteps;
                            continue;
                        }
#ifdef ENABLE_DEBUG
                        rLogger << "dF_dsigma[0] (DP) " <<  "\n" << dF_dsigma[0].transpose() << "\n" << "\n";
#endif
                    }

                    // Rounded Rankine
                    if (yieldConditionFlag(1)==ACTIVE)
                    {
                        YieldSurfaceRankine2DRoundedDerivatives(dF_dsigma[1],&(d2F_d2Sigma[1]),trialStress);
#ifdef ENABLE_DEBUG
                        rLogger << "dF_dsigma[1] (Rankine)" <<  "\n" << dF_dsigma[1].transpose() << "\n" << "\n";
                        //check second derivative
#endif
                    }


                    // ************************************************************************
                    //  residual
                    // ************************************************************************
                    residual = lastPlastStrain-rEpsilonP;

                    for (int count=0; count<numYieldSurfaces; count++)
                    {
                        if (yieldConditionFlag[count] == ACTIVE)
                        {
                            residual += deltaGamma(count)*dF_dsigma[count];
                        }
                    }

#ifdef ENABLE_DEBUG
                    rLogger << "residual " <<  "\n" << residual.transpose() << "\n" << "\n";
#endif

                    //this is just for scaling with a relative norm
                    double absResidual = residual.norm()/f_ct*e_mod;

#ifdef ENABLE_DEBUG
                    rLogger << iteration <<" residual " << absResidual << " yield condition " << yieldCondition.transpose() << "\n" << "\n";
#endif

                    // in case of PERFECT PLASTICITY [A] = hessian
                    hessian = dElastInv;

                    for (int count=0; count<numYieldSurfaces; count++)
                    {
                        if (yieldConditionFlag(count)==ACTIVE)
                        {
#ifdef ENABLE_DEBUG
                    rLogger << iteration <<" d2F_d2Sigma[" << count << "] "<< "\n" << d2F_d2Sigma[count] << "\n" << "\n";
#endif
                            hessian+=deltaGamma(count)*d2F_d2Sigma[count];
                        }
                    }
#ifdef ENABLE_DEBUG
                    rLogger << iteration <<" hessian" << "\n" << hessian << "\n" << "\n";

                    if (fabs(hessian.determinant())<toleranceDeterminant)
                    {
                        rLogger << "hessian"<< "\n" << hessian << "\n";
                        rLogger << "trialStress"<< "\n" << trialStress << "\n";
                        rLogger << "yieldConditionFlag " <<  "\n" << yieldConditionFlag.transpose() << "\n" << "\n";
                    }
#endif
                    assert(fabs(hessian.determinant())>toleranceDeterminant);

                    hessian = hessian.inverse().eval();
#ifdef ENABLE_DEBUG
                    rLogger << "determinant of hessian" << hessian.determinant() << "\n";
#endif

                    // check abs_residual and yieldCondition for active yield surfaces
                    bool convergenceFlagYieldCondition(true);
                    if (absResidual > toleranceResidual)
                        convergenceFlagYieldCondition = false;
                    else
                    {
                        for (int count=0; count<numYieldSurfaces; count++)
                        {
                            if (yieldConditionFlag(count) == INACTIVE)
                                continue;
                            if (fabs(yieldCondition(count)) > toleranceYieldSurface*f_ct)
                            {
                                convergenceFlagYieldCondition = false;
                                break;
                            }
                        }
                    }

                    if (convergenceFlagYieldCondition==true)
                    {
                        // convergence is achieved - now check of the deltaGamma is nonnegative and all other yield surfaces are valid
                        for (int count=0; count<numYieldSurfaces; count++)
                        {
                            if (yieldConditionFlag(count) == INACTIVE)
                            {
                                if (yieldCondition(count) > toleranceYieldSurface*f_ct)
                                {
                                    convergenceFlagYieldCondition = false;
                                    iteration=maxSteps;
                                    break;
                                }
                            }
                            else
                            {
                                if (deltaGamma(count)<0)
                                {
                                    convergenceFlagYieldCondition = false;
                                    iteration=maxSteps;
                                    break;
                                }
                            }
                        }

                        if (convergenceFlagYieldCondition)
                        {
                            convergedInternal = true;
#ifdef ENABLE_DEBUG
                            rLogger << "convergence after " << iteration << " iterations" << "\n" << "\n";
#endif
                        }
                    }
                    if (convergedInternal)
                    {
                        if (cutbackFactorExternal==1)
                        {
                            // compute elasto plastic d_matrix
                            // compute G_matrix
                            int curYieldFunction = 0;
                            matG.setZero(numActiveYieldFunctions,numActiveYieldFunctions);
                            vectorN.resize(numYieldSurfaces);
                            for (int count=0; count<numYieldSurfaces; count++)
                            {
                                if (yieldConditionFlag(count)==INACTIVE)
                                    continue;
                                int curYieldFunction2 = 0;
                                for (int count2=0; count2<=count; count2++)
                                {
                                    if (yieldConditionFlag(count2)==INACTIVE)
                                        continue;

                                    matG(curYieldFunction,curYieldFunction2) = (dF_dsigma[count].transpose() * hessian * dF_dsigma[count2])(0);
                                    // copy symmetric part
                                    if (count!=count2)
                                        matG(curYieldFunction2,curYieldFunction) = matG(curYieldFunction,curYieldFunction2);

                                    curYieldFunction2++;
                                }

                                // N
                                vectorN[count] = hessian * dF_dsigma[count];
#ifdef ENABLE_DEBUG
                                rLogger << "vectorN[ " << count <<"]" << "\n" << vectorN[count] <<"\n"<< "\n";
#endif
                                curYieldFunction++;
                            }

                            // solve linearized system of equations for G_inv
                            assert(fabs(matG.determinant())>toleranceDeterminant);
                            matGInv = matG.inverse();

#ifdef ENABLE_DEBUG
                        rLogger << "matG " << "\n" << matG <<"\n"<< "\n";
                        rLogger << "matGInv " << "\n" << matGInv <<"\n"<< "\n";
                        rLogger << "hessian " << "\n" << hessian <<"\n"<< "\n";
#endif
                            // compute elasto_plastic matrix
                            tmpStiffness = hessian;
                            curYieldFunction = 0;
                            for (int count=0; count<numYieldSurfaces; count++)
                            {
                                if (yieldConditionFlag(count)==INACTIVE)
                                    continue;
                                int curYieldFunction2 = 0;
                                for (int count2=0; count2<numYieldSurfaces; count2++)
                                {
                                    if (yieldConditionFlag(count2)==INACTIVE)
                                        continue;
                                    tmpStiffness-=matGInv(curYieldFunction,curYieldFunction2)*vectorN[count]*vectorN[count2].transpose();
                                    curYieldFunction2++;
                                }
                                curYieldFunction++;
                            }

                            //update new history variables
                            rdEpsilonPdEpsilon = dElastInv * (dElast - tmpStiffness);
#ifdef ENABLE_DEBUG
                        rLogger << "rdEpsilonPdEpsilon " << "\n" << rdEpsilonPdEpsilon <<"\n"<< "\n";
                        rLogger << "dElastInv " << "\n" << dElastInv <<"\n"<< "\n";
                        rLogger << "dElast " << "\n" << dElast <<"\n"<< "\n";
                        rLogger << "tmpStiffness " << "\n" << tmpStiffness <<"\n"<< "\n";
#endif

                            rStress = trialStress;

    #ifdef ENABLE_DEBUG
                            rLogger << "numberOfExternalSteps (totalCurstrainincreases) " << numberOfExternalCutbacks <<"\n";
                            rLogger << "numberOfInternalIterations (Newton iterations to find delta2Gamma and deltaStress) " << numberOfInternalIterations <<"\n";
    #endif
                        }
                        else
                        {
                            //do nothing and just increase the strain
                        }
                    }
                    else
                    {
                        //no convergence convergedInternal==false
                        if (iteration<maxSteps)
                        {
                            //compute delta gamma
                            int curYieldFunction = 0;
                            matG.setZero(numActiveYieldFunctions,numActiveYieldFunctions);
                            for (int count=0; count<numYieldSurfaces; count++)
                            {
                                if (yieldConditionFlag(count)==INACTIVE)
                                    continue;
                                int curYieldFunction2 = 0;
                                for (int count2=0; count2<=count; count2++)
                                {
                                    if (yieldConditionFlag(count2)==INACTIVE)
                                        continue;

                                    matG(curYieldFunction,curYieldFunction2) = (dF_dsigma[count].transpose()*hessian*dF_dsigma[count2])(0);
                                    // copy symmetric part
                                    if (count!=count2)
                                    {
                                        matG(curYieldFunction2,curYieldFunction) = matG(curYieldFunction,curYieldFunction2);
                                    }

                                    curYieldFunction2++;
                                }

                                curYieldFunction++;
                            }
                            assert(curYieldFunction==numActiveYieldFunctions);

                            // solve linearized system of equations for G_inv
                            assert(fabs(matG.determinant())>toleranceDeterminant);
#ifdef ENABLE_DEBUG
                            rLogger << "G and determinant" << matG.determinant() << "\n";
                            rLogger << matG << "\n";
#endif

                            matGInv = matG.inverse();

                            // compute deltaGamma
                            Eigen::Matrix<double,4,1> helpVector = hessian * residual;
                            Eigen::Matrix<double,4,1> helpVector2(0,0,0,0);
                            for (int count=0; count<numYieldSurfaces; count++)
                            {
                                if (yieldConditionFlag(count)==INACTIVE)
                                    continue;

                                helpVector2(count) = yieldCondition(count) - dF_dsigma[count].dot(helpVector);
                            }

                            curYieldFunction = 0;
                            for (int count=0; count<numYieldSurfaces; count++)
                            {
                                if (yieldConditionFlag(count)==INACTIVE)
                                    continue;
                                delta2Gamma(count) = 0.;
                                int curYieldFunction2 = 0;
                                for (int count2=0; count2<numYieldSurfaces; count2++)
                                {
                                    if (yieldConditionFlag(count2)==INACTIVE)
                                        continue;

                                    delta2Gamma(count) += matGInv(curYieldFunction,curYieldFunction2)*helpVector2(count2);

                                    curYieldFunction2++;
                                }
                                curYieldFunction++;
                            }

#ifdef ENABLE_DEBUG
                            rLogger << "delta2Gamma " << delta2Gamma.transpose() << "\n"<< "\n";
#endif
                            // ******************************************************************
                            //  compute increments for stress
                            // ******************************************************************
                            helpVector = residual;

                            for (int count=0; count<numYieldSurfaces; count++)
                            {
                                if (yieldConditionFlag(count)==INACTIVE)
                                    continue;
#ifdef ENABLE_DEBUG
                                rLogger<<"dfdsigma " << dF_dsigma[count].transpose()<<"\n"<< "\n";
#endif
                                helpVector += delta2Gamma(count)*dF_dsigma[count];
                            }

                            deltaStress =  hessian * helpVector;
#ifdef ENABLE_DEBUG
                            rLogger<<"deltaStress " << deltaStress.transpose()<<"\n"<< "\n";
#endif
                            //internal line search convergedExternal
                            double deltaCutbackFactor(1.);
                            double cutbackFactorLS(deltaCutbackFactor);
                            bool convergedLS(false);

                            //norm of the residual in the local Newton iteration, used as convergence indicator
                            // in the local iteration, the unknowns are deltaSigma and delta2Gamma
                            // whereas the residuals are the difference between (total -elastic) and (plastic strain)
                            // and the yield conditions
                            double normInit = residual.squaredNorm();
                            for (int count=0; count<numYieldSurfaces; count++)
                            {
                                if (yieldConditionFlag(count)==INACTIVE)
                                    continue;
                                normInit +=yieldCondition(count)*yieldCondition(count);
                            }
                            int numberOfLinesearchSteps(0);
                            while (!convergedLS)
                            {
                                numberOfLinesearchSteps++;
                                convergedLS = true;

                                deltaGammaLS= deltaGamma + cutbackFactorLS * delta2Gamma;
                                deltaPlasticStrainLS = cutbackFactorLS * (dElastInv * deltaStress);
                                epsilonPLS = rEpsilonP + deltaPlasticStrainLS;
                                stressLS = trialStress - cutbackFactorLS*deltaStress;

#ifdef ENABLE_DEBUG
                                rLogger << "delta2Gamma " << delta2Gamma.transpose() << "\n"<< "\n";
                                rLogger << "deltaPlasticStrainLS " << deltaPlasticStrainLS.transpose() << "\n"<< "\n";
                                rLogger << "epsilonPLS " << epsilonPLS.transpose() << "\n"<< "\n";
                                rLogger << "stressLS " << stressLS.transpose() << "\n"<< "\n";
#endif


                                // calculate yield condition and yield condition flag (active or not)
                                // Drucker Prager
                                if (yieldConditionFlag(0)==ACTIVE)
                                {
                                    yieldConditionLS(0) = YieldSurfaceDruckerPrager2D(stressLS, BETA, H_P);
                                    if (!YieldSurfaceDruckerPrager2DDerivatives(dF_dsigma[0],0,stressLS,BETA))
                                    {
                                        //no convergence, decrease line search step
                                        convergedLS =false;
                                    }
#ifdef ENABLE_DEBUG
                                    rLogger << "dF_dsigma[0] " <<  "\n" << dF_dsigma[0].transpose() << "\n" << "\n";
#endif
                                }

                                // Rounded Rankine
                                if (yieldConditionFlag(1)==ACTIVE)
                                {
                                    yieldConditionLS(1) = YieldSurfaceRankine2DRounded(stressLS, f_ct);
                                    YieldSurfaceRankine2DRoundedDerivatives(dF_dsigma[1],0,stressLS);
#ifdef ENABLE_DEBUG
                                    rLogger << "dF_dsigma[1] " <<  "\n" << dF_dsigma[1].transpose() << "\n" << "\n";
#endif
                                }

                                // residual in line search
                                residualLS = lastPlastStrain-epsilonPLS;

                                for (int count=0; count<numYieldSurfaces; count++)
                                {
                                    if (yieldConditionFlag[count] == ACTIVE)
                                    {
                                        residualLS += deltaGammaLS(count)*dF_dsigma[count];
                                    }
                                }

#ifdef ENABLE_DEBUG
                                rLogger << "residual linesearch" <<  "\n" << residualLS.transpose() << "\n" << "\n";
#endif
                                double normCurr = residualLS.squaredNorm();
                                for (int count=0; count<numYieldSurfaces; count++)
                                {
                                    if (yieldConditionFlag(count)==INACTIVE)
                                        continue;
                                    normCurr +=yieldConditionLS(count)*yieldConditionLS(count);
                                }
#ifdef ENABLE_DEBUG
                                rLogger << "normInit " << normInit << "normCurr " << normCurr<<"\n"<< "\n";
#endif

                                // check of relative norm of residual is sufficiently decreasing or if it is almost zero
                                if ((normCurr/normInit)*(normCurr/normInit)>1-0.5*cutbackFactorLS && normCurr>toleranceResidual*toleranceResidual)
                                    convergedLS = false;
                                if (cutbackFactorLS<=minCutbackFactorLS)
                                    convergedLS = true;

                                if (!convergedLS)
                                {
                                    deltaCutbackFactor*=0.5;
                                    cutbackFactorLS-=deltaCutbackFactor;
                                }
                            }
                            trialStress = stressLS;
                            deltaGamma = deltaGammaLS;
                            rEpsilonP = epsilonPLS;

#ifdef ENABLE_DEBUG
                            if (yieldConditionFlag(0)==ACTIVE)
                                rLogger << "dF_dsigma[0] at end of line search " <<  "\n" << dF_dsigma[0].transpose() << "\n" << "\n";
                            if (yieldConditionFlag(1)==ACTIVE)
                                rLogger << "dF_dsigma[1] at end of line search " <<  "\n" << dF_dsigma[1].transpose() << "\n" << "\n";
                            rLogger << "numberOfLinesearchSteps " << numberOfLinesearchSteps << "\n";
#endif
                        }
                    }
                } // end of loop
#ifdef ENABLE_DEBUG
                if (convergedInternal==false)
                {
                    rLogger << "state with fixed yield conditions did not converge. norm Residual " << residual.squaredNorm() << "\n";
                    for (int count=0; count<2; count++)
                    {
                        if (yieldConditionFlag(count)==ACTIVE)
                        {
                            rLogger << "     yield condition "<< count+1 << " " << yieldCondition(count) << " ("<<toleranceYieldSurface*f_ct << ")"<<"\n";
                            rLogger << "     deltaGamma "<< count+1 << " " << deltaGamma(count) << "\n";
                        }
                        else
                            rLogger << "     yield condition "<< count+1 << " " << yieldCondition(count) << "\n";

                    }
                }
#endif
            }
        }
        catch (NuTo::MechanicsException& e)
        {
            e.AddMessage("[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] Error performing return mapping procedure.");
            throw e;
        }
        catch (...)
        {
               throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] Error performing return mapping procedure.");
        }

        if (convergedInternal)
        {
#ifdef ENABLE_DEBUG
            rLogger << "numberOfInternalIterations " << numberOfInternalIterations -  prevNumberOfInternalIterations<< "(" << numberOfInternalIterations << ")" << "\n";
            rLogger << "convergence for external cutback factor" << "\n";
            prevNumberOfInternalIterations = numberOfInternalIterations;
#endif

            //update equivalente plastic strain
            deltaPlasticStrain = rEpsilonP - lastPlastStrain;
            rDeltaEqPlasticStrain += sqrt(deltaPlasticStrain[0]*deltaPlasticStrain(0)+deltaPlasticStrain(1)*deltaPlasticStrain(1)+
                            0.5 *(deltaPlasticStrain(2)*deltaPlasticStrain(2))+deltaPlasticStrain(3)*deltaPlasticStrain(3));

//            double tmp(sqrt(rEpsilonP(0)*rEpsilonP(0) + rEpsilonP(1)*rEpsilonP(1)+ 0.5*rEpsilonP(2)*rEpsilonP(2)+ rEpsilonP(3)*rEpsilonP(3)));
//            rLogger << "\n" << "rDeltaEqPlasticStrain " << rDeltaEqPlasticStrain << "Norm of epsilonp " << tmp << "delta between " << rDeltaEqPlasticStrain - tmp << "\n"<< "\n";
//            rLogger << "\n" << "rEpsilonP " << rEpsilonP << "Norm of epsilonp " << "\n";

            if (cutbackFactorExternal==1.)
            {
                convergedExternal=true;
            }
            else
            {
                lastPlastStrain = rEpsilonP;
                lastDeltaEqPlasticStrain = rDeltaEqPlasticStrain;

                if (numberOfInternalIterations<10)
                    deltaCutbackFactorExternal*=1.5;
                cutbackFactorExternal +=deltaCutbackFactorExternal;
                if (cutbackFactorExternal>1)
                {
                    deltaCutbackFactorExternal -= cutbackFactorExternal -1.;
                    cutbackFactorExternal = 1.;
                }
            }
        }
        else
        {
            // no convergence in return mapping
            deltaCutbackFactorExternal*=0.5;
            cutbackFactorExternal-=deltaCutbackFactorExternal;
            //#ifdef ENABLE_DEBUG
                        rLogger << "decrease external cutback factor to " << deltaCutbackFactorExternal << "\n";
            //#endif
        }

    }
    if (cutbackFactorExternal<=minCutbackFactor)
    {
        rLogger << "[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] No convergence can be obtained in the return mapping procedure." << "n";
        return Error::NO_CONVERGENCE;
    }
*/
    return Error::SUCCESSFUL;
}

//! @brief ... performs the return mapping procedure for the plasticity model
//! @param rStrain              ... current total strain
//! @param rPrevPlasticStrain   ... previous plastic strain (history variable)
//! @param rPrevTotalStrain     ... previous total strain (history variable)
//! @param rPrevEqPlasticStrain ... previous equiavalente plastic strain (history variable)
//! @param rEpsilonP            ... new plastic strain after return mapping
//! @param rEqPlasticStrain     ... new equivalente olastic strain after return mapping
//! @param rdSigmadEpsilon      ... new derivative of current stress with respect to the total strain
//! @param rdEpsilonPdEpsilon   ... new derivative of current plastic strain with respect to the total strain
NuTo::Error::eError NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping3D(
        const EngineeringStrain3D& rStrain,
        const double rPrevPlasticStrain[6],
        const EngineeringStrain3D& rPrevTotalStrain,
        Eigen::Matrix<double,6,1>& rStress,
        Eigen::Matrix<double,6,1>& rEpsilonP,
        double& rDeltaEqPlasticStrain,
        Eigen::Matrix<double,6,6>& rdSigmadEpsilon,
        Eigen::Matrix<double,6,6>& rdEpsilonPdEpsilon,
        NuTo::Logger& rLogger)const
{

	double e_mod = mE; //modify that one in the case of random fields
	double   nu  = mNu;
	double f_ct  = mTensileStrength;
	double f_c1  = mCompressiveStrength;
	double f_c2  = mBiaxialCompressiveStrength;

	assert(f_c2>f_c1);
	assert(f_c1>0);
	assert(f_c2>0);
	assert(e_mod>0);

	// ******************************************************************
	// *    F_BETA:    required by DRUCKER-PRAGER yield surface               *
	// ******************************************************************
	double BETA = sqrt3*(f_c2-f_c1) / (2*f_c2-f_c1);
	double H_P  = f_c2*f_c1 / (sqrt3*(2*f_c2-f_c1));

	//! @brief strain currently solved for the plastic strains, in general equal to  rStrain, but if applied in steps, it's smaller
	Eigen::Matrix<double,6,1> curTotalStrain;
	//! @brief previous plastic strain, either from previous equilibrium (static data) or if applied in steps, previous converged state
	Eigen::Matrix<double,6,1> lastPlastStrain;
	//! @brief plastic strain in the line search
	Eigen::Matrix<double,6,1> epsilonPLS;
	//! @brief previous eq plastic strain, either from previous equilibrium (static data) or if applied in steps, previous converged state
	double lastDeltaEqPlasticStrain;
	//! @brief residual in the return mapping procedure
	Eigen::Matrix<double,6,1> residual;
	//! @brief residual in the return mapping procedure within linesearch
	Eigen::Matrix<double,6,1> residualLS;
	//! @brief full stress increment within one iteration of the return mapping, might be applied in steps in the succeeding line search
	Eigen::Matrix<double,6,1> deltaStress;
	//! @brief total strain increment between strain from previous static data and new total strain
	Eigen::Matrix<double,6,1> deltaStrain;
	//! @brief plastic strain increment within the linesearch
	Eigen::Matrix<double,6,1> deltaPlasticStrainLS;
	//! @brief plastic strain increment between to load steps (internal load splitting)
	Eigen::Matrix<double,6,1> deltaPlasticStrain;
	//! @brief trial stress of the first iteration
	Eigen::Matrix<double,6,1> initTrialStress;
	//! @brief trial stress in the line search
	Eigen::Matrix<double,6,1> stressLS;
	//! @brief trial stress
	Eigen::Matrix<double,6,1> trialStress;
	//! @brief elastic strain
	Eigen::Matrix<double,6,1> elasticStrain;
	//! @brief elastic strain in line search
	Eigen::Matrix<double,6,1> elasticStrainLS;
	//! @brief plastic multiplier
	Eigen::Matrix<double,2,1> deltaGamma;
	//! @brief plastic multiplier
	Eigen::Matrix<double,2,1> deltaGammaLS;
	//! @brief increment of plastic multiplier in return mapping procedure
	Eigen::Matrix<double,2,1> delta2Gamma;
	//! @brief yield condition
	Eigen::Matrix<double,2,1> yieldCondition;
	//! @brief yield condition in line search
	Eigen::Matrix<double,2,1> yieldConditionLS(0,0);
	//! @brief yield condition at the first iteration
	Eigen::Matrix<double,2,1> initYieldCondition;
	//! @brief flag that indicates if a yield function is active or not
	Eigen::Matrix<bool,2,1> yieldConditionFlag;
	//! @brief (dF/dsigma)T * Hessian * dF/dsigma
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matG;
	//! @brief ((dF/dsigma)T * Hessian * dF/dsigma )^-1
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matGInv;
	//! @brief algorithmic modulus DElastInv+delta_gamm*d2F/d2Sigma
	Eigen::Matrix<double,6,6> hessian;
	//! @brief first derivatives of the yield functions
	std::vector<Eigen::Matrix<double,6,1> > dF_dsigma;
	//! @brief second derivatives of the yield functions
	std::vector<Eigen::Matrix<double,6,6> > d2F_d2Sigma;
	//! @brief algorithmic modulus * dF_dsigma
	std::vector<Eigen::Matrix<double,6,1> > vectorN;
	//! @brief number of active yield functions
	int numActiveYieldFunctions;

	// for the application of strains in steps, calculate the total strain increment to be applied
	deltaStrain(0) = rStrain.mEngineeringStrain[0]-rPrevTotalStrain.mEngineeringStrain[0];
	deltaStrain(1) = rStrain.mEngineeringStrain[1]-rPrevTotalStrain.mEngineeringStrain[1];
	deltaStrain(2) = rStrain.mEngineeringStrain[2]-rPrevTotalStrain.mEngineeringStrain[2];
	deltaStrain(3) = rStrain.mEngineeringStrain[3]-rPrevTotalStrain.mEngineeringStrain[4];
	deltaStrain(4) = rStrain.mEngineeringStrain[4]-rPrevTotalStrain.mEngineeringStrain[5];
	deltaStrain(5) = rStrain.mEngineeringStrain[5]-rPrevTotalStrain.mEngineeringStrain[6];
#ifdef ENABLE_DEBUG
	rLogger << "\n" << "deltaStrain" << deltaStrain.transpose() << "\n" << "\n";
#endif

	// initialize last plastic strain and last converged stress
	lastPlastStrain << rPrevPlasticStrain[0], rPrevPlasticStrain[1] ,rPrevPlasticStrain[2] ,
			           rPrevPlasticStrain[3], rPrevPlasticStrain[4] ,rPrevPlasticStrain[5];
#ifdef ENABLE_DEBUG
		rLogger << "\n" << "lastPlastStrain" << lastPlastStrain.transpose() << "\n" << "\n";
#endif
	rDeltaEqPlasticStrain = 0;
	// *****************************************************************
	//                   elastic matrix generation                     *
	// *****************************************************************
	//! @brief elastic stiffness
	Eigen::Matrix<double,6,6> dElast;
	//! @brief inverse elastic stiffness
	Eigen::Matrix<double,6,6> dElastInv;

	{
		double factor = e_mod/((1.+nu)*(1.-2.*nu));
		double oneminusnufactor = (1-nu)*factor;
		double nufactor = nu*factor;

		dElast <<  oneminusnufactor , nufactor         , nufactor , 0.             , 0.              , 0.         ,
				   nufactor         , oneminusnufactor , nufactor , 0.             , 0.              , 0.         ,
				   nufactor         , nufactor , oneminusnufactor , 0.             , 0.              , 0.         ,
				   0.               , 0.               , 0.       ,(0.5-nu)*factor , 0.              , 0.         ,
				   0.               , 0.               , 0.       , 0.             , (0.5-nu)*factor , 0.,
				   0.               , 0.               , 0.       , 0.             , 0.              , (0.5-nu)*factor;

		factor = 1./e_mod;
		double minusnufactor = -nu*factor;
		dElastInv   << factor        , minusnufactor , minusnufactor,0.                , 0.               , 0.              ,
					   minusnufactor , factor        , minusnufactor,0.                , 0.               , 0.              ,
					   minusnufactor , minusnufactor , factor       ,0.                , 0.               , 0.              ,
					   0.            , 0.            , 0.           , 2.*factor*(1.+nu), 0.               , 0.              ,
					   0.            , 0.            , 0.           , 0.               , 2.*factor*(1.+nu), 0.              ,
					   0.            , 0.            , 0.           , 0.               , 0.               ,2.*factor*(1.+nu);

	}
	//! @brief delta load factor for the previous iteration
	double deltaCutbackFactorExternal(1.);

	//! @brief current load factor (between 0 and 1) to apply the total strain increment in steps
	double cutbackFactorExternal(deltaCutbackFactorExternal);

	//! @brief flag to determine if the iteration is finished (converged at  cutbackFactorExternal=1)
	bool convergedExternal(false);

	int numberOfExternalCutbacks(0);
	int numberOfInternalIterations(0);
#ifdef ENABLE_DEBUG
	int prevNumberOfInternalIterations(0);
#endif
	lastDeltaEqPlasticStrain = 0.;
	while (cutbackFactorExternal>minCutbackFactor && !convergedExternal)
	{
		numberOfExternalCutbacks++;

		curTotalStrain(0) = rPrevTotalStrain.mEngineeringStrain[0]+cutbackFactorExternal*deltaStrain(0);
		curTotalStrain(1) = rPrevTotalStrain.mEngineeringStrain[1]+cutbackFactorExternal*deltaStrain(1);
		curTotalStrain(2) = rPrevTotalStrain.mEngineeringStrain[2]+cutbackFactorExternal*deltaStrain(2);
		curTotalStrain(3) = rPrevTotalStrain.mEngineeringStrain[3]+cutbackFactorExternal*deltaStrain(3);
		curTotalStrain(4) = rPrevTotalStrain.mEngineeringStrain[4]+cutbackFactorExternal*deltaStrain(4);
		curTotalStrain(5) = rPrevTotalStrain.mEngineeringStrain[5]+cutbackFactorExternal*deltaStrain(5);

#ifdef ENABLE_DEBUG
		rLogger << "\n" << "curTotalStrain " << curTotalStrain.transpose() << "\n";
//		rLogger << "\n" << "rPrevTotalStrain" << rPrevTotalStrain << "\n";
		rLogger << "\n" << "cutbackFactorExternal " << cutbackFactorExternal << "\n";
#endif

		// checks the convergence of the Newton iteration for a prescribed current strain
		bool convergedInternal(false);
		try
		{
			//resize yield condition vector
			int numYieldSurfaces=2;
			yieldCondition.setZero(numYieldSurfaces);

			//elastic strain, stress and d_matrix
			elasticStrain = curTotalStrain - lastPlastStrain;

			//TODO just use the upper part for new EigenVersion 3.0
			trialStress = dElast.selfadjointView<Eigen::Upper>()*elasticStrain;
			initTrialStress = dElast*elasticStrain;
#ifdef ENABLE_DEBUG
			rLogger << "initTrialStress " << "\n" << initTrialStress.transpose() << "\n" << "\n";
#endif

			//calculate yield condition
			//Drucker Prager
			bool errorDerivative(false);
			initYieldCondition(0) = YieldSurfaceDruckerPrager3D(initTrialStress, BETA, H_P, 0 ,0, errorDerivative);

			//rounded Rankine (0,0 no first and second derivative
			initYieldCondition(1) = YieldSurfaceRankine3DRounded(initTrialStress, f_ct, 0, 0);

#ifdef ENABLE_DEBUG
			rLogger << "initYieldCondition " << "\n" << initYieldCondition(0) << " " << initYieldCondition(1) << "\n";
#endif


			if (initYieldCondition(0)<-toleranceYieldSurface*f_ct && initYieldCondition(1)<-toleranceYieldSurface*f_ct)
			{
				// *************************************************
				// *  thus we have elastic -------------> elastic  *
				// *************************************************
#ifdef ENABLE_DEBUG
				rLogger << "linear elastic step" << "\n" << "\n";
#endif
				convergedInternal = true;
				rEpsilonP =  lastPlastStrain;
				rDeltaEqPlasticStrain = 0;
				trialStress = initTrialStress;
				if (cutbackFactorExternal==1)
				{
					rdEpsilonPdEpsilon.setZero(6,6);
					rStress = trialStress;
					return Error::SUCCESSFUL;
				}
			}
			else
			{
#ifdef ENABLE_DEBUG
				rLogger << "plastic step" << "\n" << "\n";
#endif
			}

			// perform return mapping
			dF_dsigma.resize(numYieldSurfaces);
			d2F_d2Sigma.resize(numYieldSurfaces);

			// initialize plastic multiplier
			deltaGamma.setZero(numYieldSurfaces);
			delta2Gamma.setZero(numYieldSurfaces);
			yieldConditionFlag.setZero(numYieldSurfaces);

			for (int fixedYieldConditions=0; fixedYieldConditions<3 && convergedInternal==false; fixedYieldConditions++)
			{
				switch(fixedYieldConditions)
				{
				case 0:
#ifdef ENABLE_DEBUG
					rLogger<< "pure Rankine" << "\n";
#endif
					deltaGamma(1) = 0;
					if (initYieldCondition(1)<-toleranceYieldSurface*f_ct)
						continue;
					yieldConditionFlag(0) = INACTIVE;
					yieldConditionFlag(1) = ACTIVE;
					numActiveYieldFunctions = 1;
					break;
				case 1:
#ifdef ENABLE_DEBUG
					rLogger<< "combined" << "\n";
#endif
					if (initYieldCondition(0)<-toleranceYieldSurface*f_ct || initYieldCondition(1)<-toleranceYieldSurface*f_ct)
						continue;
					yieldConditionFlag(0) = ACTIVE;
					yieldConditionFlag(1) = ACTIVE;
					deltaGamma(0) = 0;
					deltaGamma(1) = 0;
					numActiveYieldFunctions = 2;
					break;
				case 2:
#ifdef ENABLE_DEBUG
					rLogger<< "pure DP" << "\n";
#endif
					if (initYieldCondition(0)<-toleranceYieldSurface*f_ct)
						continue;
					yieldConditionFlag(0) = ACTIVE;
					yieldConditionFlag(1) = INACTIVE;
					deltaGamma(0) = 0;
					numActiveYieldFunctions = 1;
					break;
				default:
					throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] programming error - should not happen.");
				}

				for (int iteration = 0; iteration < maxSteps && convergedInternal==false; iteration++)
				{
					numberOfInternalIterations++;
					if (iteration==0)
					{
						rEpsilonP = lastPlastStrain;
						rDeltaEqPlasticStrain = lastDeltaEqPlasticStrain;
						yieldCondition= initYieldCondition;
						trialStress = initTrialStress;
					}
					else
					{
						//trial stress is the last stress state from the previous line search
						if (yieldConditionFlag(0)==ACTIVE)
							yieldCondition(0) = yieldConditionLS(0);
						else
							yieldCondition(0) = YieldSurfaceDruckerPrager3D(trialStress, BETA, H_P, 0, 0, errorDerivative);

						if (yieldConditionFlag(1)==ACTIVE)
							yieldCondition(1) = yieldConditionLS(1);
						else
							yieldCondition(1) = YieldSurfaceRankine3DRounded(trialStress, f_ct, 0, 0);

					}
#ifdef ENABLE_DEBUG
					rLogger << "trialStress " <<  "\n" << trialStress.transpose() << "\n" << "\n";
					rLogger << "yieldCondition " <<  "\n" << yieldCondition.transpose() << "\n" << "\n";
#endif

					// DP
					if (yieldConditionFlag(0)==ACTIVE)
					{
						YieldSurfaceDruckerPrager3D(trialStress, BETA, H_P,&(dF_dsigma[0]),&(d2F_d2Sigma[0]), errorDerivative);
                        if (errorDerivative)
						{
                            //no convergence, decrease line search step
                            iteration = maxSteps;
                            continue;
                        }
#ifdef ENABLE_DEBUG
						rLogger << "dF_dsigma[0] (DP) " <<  "\n" << dF_dsigma[0].transpose() << "\n" << "\n";
#endif
					}

					// Rounded Rankine
					if (yieldConditionFlag(1)==ACTIVE)
					{
						YieldSurfaceRankine3DRounded(trialStress,f_ct, &(dF_dsigma[1]),&(d2F_d2Sigma[1]));
#ifdef ENABLE_DEBUG
						rLogger << "dF_dsigma[1] (Rankine)" <<  "\n" << dF_dsigma[1].transpose() << "\n" << "\n";
						//check the calculation
						double delta(1e-8);
						//if (rd2F_d2Sigma!=0)
						{
						//rLogger << "\n" << "check yield surface and derivatives" << "\n";
							//rLogger << "sigmas " << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
							Eigen::Matrix<double,6,1>stress(trialStress);
							Eigen::Matrix<double,6,1> dF_dSigma1, dF_dSigma2,dF_dSigmaCDF;
							Eigen::Matrix<double,6,6> d2F_d2SigmaExact,d2F_d2SigmaCDF;
							double f1 = YieldSurfaceRankine3DRounded(rStress,f_ct,&dF_dSigma1,&d2F_d2SigmaExact);
							for (int count=0; count<6; count++)
							{
								stress(count)+= delta;
								double f2 = YieldSurfaceRankine3DRounded(stress,f_ct,&dF_dSigma2,0);
								dF_dSigmaCDF(count) = (f2-f1)/delta;

								d2F_d2SigmaCDF.row(count) = (dF_dSigma2-dF_dSigma1)/delta;

								stress(count)-=delta;
							}

							if ((dF_dSigmaCDF-dF_dSigma1).array().abs().maxCoeff()>1e-1)
							{
								//std::cout << "sigmas principal" << principal[0]<< " " << principal[1] <<  " " << principal[2] << "\n";
								std::cout << "sigmas " << stress(0) << " " << stress(1) <<  " " << stress(2) <<  " " <<stress(3) << "\n";
								std::cout << "error first derivative " << (dF_dSigmaCDF-dF_dSigma1).array().abs().maxCoeff() << "\n";

								std::cout<< "rdF_dSigma " << "\n" << dF_dSigma1 << "\n"<< "\n";
								std::cout<< "rdF_dSigmaCDF " << "\n" << dF_dSigmaCDF << "\n"<< "\n";
								//throw MechanicsException("[NuTo::Mechanics::GradientDamagePlasticityEngineeringStress] Error calculating first derivative of yield function.");
							}

							// fabs is checked, since of the type of the yield surface changes, the second derivatives are likely to change as well
							//if (fabs(principal[0])>delta && fabs(principal[1])>delta && fabs(principal[2])>delta)
							{
								if ((d2F_d2SigmaCDF-(d2F_d2SigmaExact)).array().abs().maxCoeff()>1e-1)
								{
									//std::cout << "sigmas principal" << principal[0]<< " " << principal[1] <<  " " << principal[2] << "\n";
									//std::cout << "sigmas " << rStress(0) << " " << rStress(1) <<  " " << rStress(2) <<  " "<< rStress(3) <<  " "<< rStress(4) <<  " "<< rStress(5) << "\n";
									std::cout << "error second derivatives " << (d2F_d2SigmaCDF-(d2F_d2SigmaExact)).array().abs().maxCoeff() << "\n";

									std::cout<< "rd2F_d2SigmaCDF " << "\n" << d2F_d2SigmaCDF << "\n"<< "\n";
									std::cout<< "rd2F_d2SigmaExact " << "\n" << d2F_d2SigmaExact << "\n"<< "\n";
									//throw MechanicsException("[NuTo::Mechanics::GradientDamagePlasticityEngineeringStress] Error calculating second derivative of yield function.");
								}
							}
							/*else
							{
								std::cout << "at least one principal stress close to zero, stiffness matrix is not reliable" << "\n";
							}
							*/
						}
#endif
					}


					// ************************************************************************
					//  residual
					// ************************************************************************
					residual = lastPlastStrain-rEpsilonP;

					for (int count=0; count<numYieldSurfaces; count++)
					{
						if (yieldConditionFlag[count] == ACTIVE)
						{
							residual += deltaGamma(count)*dF_dsigma[count];
						}
					}

#ifdef ENABLE_DEBUG
					rLogger << "residual " <<  "\n" << residual.transpose() << "\n" << "\n";
#endif

					//this is just for scaling with a relative norm
					double absResidual = residual.norm()/f_ct*e_mod;

#ifdef ENABLE_DEBUG
					rLogger << iteration <<" residual " << absResidual << " yield condition " << yieldCondition.transpose() << "\n" << "\n";
#endif

					// in case of PERFECT PLASTICITY [A] = hessian
					hessian = dElastInv;

					for (int count=0; count<numYieldSurfaces; count++)
					{
						if (yieldConditionFlag(count)==ACTIVE)
						{
#ifdef ENABLE_DEBUG
					rLogger << iteration <<" d2F_d2Sigma[" << count << "] "<< "\n" << d2F_d2Sigma[count] << "\n" << "\n";
#endif
							hessian+=deltaGamma(count)*d2F_d2Sigma[count];
						}
					}
#ifdef ENABLE_DEBUG
					rLogger << iteration <<" hessian" << "\n" << hessian << "\n" << "\n";

					if (fabs(hessian.determinant())<toleranceDeterminant)
					{
						rLogger << "hessian"<< "\n" << hessian << "\n";
						rLogger << "trialStress"<< "\n" << trialStress << "\n";
						rLogger << "yieldConditionFlag " <<  "\n" << yieldConditionFlag.transpose() << "\n" << "\n";
					}
#endif
					assert(fabs(hessian.determinant())>toleranceDeterminant);

					hessian = hessian.inverse().eval();
#ifdef ENABLE_DEBUG
					rLogger << "determinant of hessian" << hessian.determinant() << "\n";
#endif

					// check abs_residual and yieldCondition for active yield surfaces
					bool convergenceFlagYieldCondition(true);
					if (absResidual > toleranceResidual)
						convergenceFlagYieldCondition = false;
					else
					{
						for (int count=0; count<numYieldSurfaces; count++)
						{
							if (yieldConditionFlag(count) == INACTIVE)
								continue;
							if (fabs(yieldCondition(count)) > toleranceYieldSurface*f_ct)
							{
								convergenceFlagYieldCondition = false;
								break;
							}
						}
					}

					if (convergenceFlagYieldCondition==true)
					{
						// convergence is achieved - now check of the deltaGamma is nonnegative and all other yield surfaces are valid
						for (int count=0; count<numYieldSurfaces; count++)
						{
							if (yieldConditionFlag(count) == INACTIVE)
							{
								if (yieldCondition(count) > toleranceYieldSurface*f_ct)
								{
									convergenceFlagYieldCondition = false;
									iteration=maxSteps;
									break;
								}
							}
							else
							{
								if (deltaGamma(count)<0)
								{
									convergenceFlagYieldCondition = false;
									iteration=maxSteps;
									break;
								}
							}
						}

						if (convergenceFlagYieldCondition)
						{
							convergedInternal = true;
#ifdef ENABLE_DEBUG
							rLogger << "convergence after " << iteration << " iterations" << "\n" << "\n";
#endif
						}
					}
					if (convergedInternal)
					{
						if (cutbackFactorExternal==1)
						{
							// compute elasto plastic d_matrix
							// compute G_matrix
							int curYieldFunction = 0;
							matG.setZero(numActiveYieldFunctions,numActiveYieldFunctions);
							vectorN.resize(numYieldSurfaces);
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
								int curYieldFunction2 = 0;
								for (int count2=0; count2<=count; count2++)
								{
									if (yieldConditionFlag(count2)==INACTIVE)
										continue;

									matG(curYieldFunction,curYieldFunction2) = (dF_dsigma[count].transpose() * hessian * dF_dsigma[count2])(0);
									// copy symmetric part
									if (count!=count2)
										matG(curYieldFunction2,curYieldFunction) = matG(curYieldFunction,curYieldFunction2);

									curYieldFunction2++;
								}

								// N
								vectorN[count] = hessian * dF_dsigma[count];
#ifdef ENABLE_DEBUG
								rLogger << "vectorN[ " << count <<"]" << "\n" << vectorN[count] <<"\n"<< "\n";
#endif
								curYieldFunction++;
							}

							// solve linearized system of equations for G_inv
							assert(fabs(matG.determinant())>toleranceDeterminant);
							matGInv = matG.inverse();

#ifdef ENABLE_DEBUG
						rLogger << "matG " << "\n" << matG <<"\n"<< "\n";
						rLogger << "matGInv " << "\n" << matGInv <<"\n"<< "\n";
						rLogger << "hessian " << "\n" << hessian <<"\n"<< "\n";
#endif
							// compute elasto_plastic matrix
						    rdSigmadEpsilon = hessian;
							curYieldFunction = 0;
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
								int curYieldFunction2 = 0;
								for (int count2=0; count2<numYieldSurfaces; count2++)
								{
									if (yieldConditionFlag(count2)==INACTIVE)
										continue;
									rdSigmadEpsilon-=matGInv(curYieldFunction,curYieldFunction2)*vectorN[count]*vectorN[count2].transpose();
									curYieldFunction2++;
								}
								curYieldFunction++;
							}

							//update new history variables
							rdEpsilonPdEpsilon = dElastInv * (dElast - rdSigmadEpsilon);
#ifdef ENABLE_DEBUG
						rLogger << "rdEpsilonPdEpsilon " << "\n" << rdEpsilonPdEpsilon <<"\n"<< "\n";
						rLogger << "dElastInv " << "\n" << dElastInv <<"\n"<< "\n";
						rLogger << "dElast " << "\n" << dElast <<"\n"<< "\n";
						rLogger << "rdSigmadEpsilon " << "\n" << rdSigmadEpsilon <<"\n"<< "\n";
#endif

							rStress = trialStress;

	#ifdef ENABLE_DEBUG
							rLogger << "numberOfExternalSteps (totalCurstrainincreases) " << numberOfExternalCutbacks <<"\n";
							rLogger << "numberOfInternalIterations (Newton iterations to find delta2Gamma and deltaStress) " << numberOfInternalIterations <<"\n";
	#endif
						}
						else
						{
							//do nothing and just increase the strain
						}
					}
					else
					{
						//no convergence convergedInternal==false
						if (iteration<maxSteps)
						{
							//compute delta gamma
							int curYieldFunction = 0;
							matG.setZero(numActiveYieldFunctions,numActiveYieldFunctions);
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
								int curYieldFunction2 = 0;
								for (int count2=0; count2<=count; count2++)
								{
									if (yieldConditionFlag(count2)==INACTIVE)
										continue;

									matG(curYieldFunction,curYieldFunction2) = (dF_dsigma[count].transpose()*hessian*dF_dsigma[count2])(0);
									// copy symmetric part
									if (count!=count2)
									{
										matG(curYieldFunction2,curYieldFunction) = matG(curYieldFunction,curYieldFunction2);
									}

									curYieldFunction2++;
								}

								curYieldFunction++;
							}
							assert(curYieldFunction==numActiveYieldFunctions);

							// solve linearized system of equations for G_inv
							assert(fabs(matG.determinant())>toleranceDeterminant);
#ifdef ENABLE_DEBUG
							rLogger << "G and determinant" << matG.determinant() << "\n";
							rLogger << matG << "\n";
#endif

							matGInv = matG.inverse();

							// compute deltaGamma
							Eigen::Matrix<double,6,1> helpVector = hessian * residual;
							Eigen::Matrix<double,6,1> helpVector2; helpVector2.Zero();
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;

								helpVector2(count) = yieldCondition(count) - dF_dsigma[count].dot(helpVector);
							}

							curYieldFunction = 0;
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
								delta2Gamma(count) = 0.;
								int curYieldFunction2 = 0;
								for (int count2=0; count2<numYieldSurfaces; count2++)
								{
									if (yieldConditionFlag(count2)==INACTIVE)
										continue;

									delta2Gamma(count) += matGInv(curYieldFunction,curYieldFunction2)*helpVector2(count2);

									curYieldFunction2++;
								}
								curYieldFunction++;
							}

#ifdef ENABLE_DEBUG
							rLogger << "delta2Gamma " << delta2Gamma.transpose() << "\n"<< "\n";
#endif
							// ******************************************************************
							//  compute increments for stress
							// ******************************************************************
							helpVector = residual;

							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
#ifdef ENABLE_DEBUG
								rLogger<<"dfdsigma " << dF_dsigma[count].transpose()<<"\n"<< "\n";
#endif
								helpVector += delta2Gamma(count)*dF_dsigma[count];
							}

							deltaStress =  hessian * helpVector;
#ifdef ENABLE_DEBUG
							rLogger<<"deltaStress " << deltaStress.transpose()<<"\n"<< "\n";
#endif
							//internal line search convergedExternal
							double deltaCutbackFactor(1.);
							double cutbackFactorLS(deltaCutbackFactor);
							bool convergedLS(false);

							//norm of the residual in the local Newton iteration, used as convergence indicator
							// in the local iteration, the unknowns are deltaSigma and delta2Gamma
							// whereas the residuals are the difference between (total -elastic) and (plastic strain)
							// and the yield conditions
							double normInit = residual.squaredNorm();
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
								normInit +=yieldCondition(count)*yieldCondition(count);
							}
							int numberOfLinesearchSteps(0);
							while (!convergedLS)
							{
								numberOfLinesearchSteps++;
								convergedLS = true;

								deltaGammaLS= deltaGamma + cutbackFactorLS * delta2Gamma;
								deltaPlasticStrainLS = cutbackFactorLS * (dElastInv * deltaStress);
								epsilonPLS = rEpsilonP + deltaPlasticStrainLS;
								stressLS = trialStress - cutbackFactorLS*deltaStress;

#ifdef ENABLE_DEBUG
								rLogger << "delta2Gamma " << delta2Gamma.transpose() << "\n"<< "\n";
								rLogger << "deltaPlasticStrainLS " << deltaPlasticStrainLS.transpose() << "\n"<< "\n";
								rLogger << "epsilonPLS " << epsilonPLS.transpose() << "\n"<< "\n";
								rLogger << "stressLS " << stressLS.transpose() << "\n"<< "\n";
#endif


								// calculate yield condition and yield condition flag (active or not)
								// Drucker Prager
								if (yieldConditionFlag(0)==ACTIVE)
								{
									yieldConditionLS(0) = YieldSurfaceDruckerPrager3D(stressLS, BETA, H_P,&(dF_dsigma[0]),0,errorDerivative);
									if (errorDerivative)
									{
										//no convergence, decrease line search step
										convergedLS =false;
									}
#ifdef ENABLE_DEBUG
									rLogger << "dF_dsigma[0] " <<  "\n" << dF_dsigma[0].transpose() << "\n" << "\n";
#endif
								}

								// Rounded Rankine
								if (yieldConditionFlag(1)==ACTIVE)
								{
									yieldConditionLS(1) = YieldSurfaceRankine3DRounded(stressLS, f_ct,&(dF_dsigma[1]),0);
#ifdef ENABLE_DEBUG
									rLogger << "dF_dsigma[1] " <<  "\n" << dF_dsigma[1].transpose() << "\n" << "\n";
#endif
								}

								// residual in line search
								residualLS = lastPlastStrain-epsilonPLS;

								for (int count=0; count<numYieldSurfaces; count++)
								{
									if (yieldConditionFlag[count] == ACTIVE)
									{
										residualLS += deltaGammaLS(count)*dF_dsigma[count];
									}
								}

#ifdef ENABLE_DEBUG
								rLogger << "residual linesearch" <<  "\n" << residualLS.transpose() << "\n" << "\n";
#endif
								double normCurr = residualLS.squaredNorm();
								for (int count=0; count<numYieldSurfaces; count++)
								{
									if (yieldConditionFlag(count)==INACTIVE)
										continue;
									normCurr +=yieldConditionLS(count)*yieldConditionLS(count);
								}
#ifdef ENABLE_DEBUG
								rLogger << "normInit " << normInit << "normCurr " << normCurr<<"\n"<< "\n";
#endif

								// check of relative norm of residual is sufficiently decreasing or if it is almost zero
								if ((normCurr/normInit)*(normCurr/normInit)>1-0.5*cutbackFactorLS && normCurr>toleranceResidual*toleranceResidual)
									convergedLS = false;
								if (cutbackFactorLS<=minCutbackFactorLS)
									convergedLS = true;

								if (!convergedLS)
								{
									deltaCutbackFactor*=0.5;
									cutbackFactorLS-=deltaCutbackFactor;
								}
							}
							trialStress = stressLS;
							deltaGamma = deltaGammaLS;
							rEpsilonP = epsilonPLS;

#ifdef ENABLE_DEBUG
							if (yieldConditionFlag(0)==ACTIVE)
								rLogger << "dF_dsigma[0] at end of line search " <<  "\n" << dF_dsigma[0].transpose() << "\n" << "\n";
							if (yieldConditionFlag(1)==ACTIVE)
								rLogger << "dF_dsigma[1] at end of line search " <<  "\n" << dF_dsigma[1].transpose() << "\n" << "\n";
							rLogger << "numberOfLinesearchSteps " << numberOfLinesearchSteps << "\n";
#endif
						}
					}
				} // end of loop
#ifdef ENABLE_DEBUG
				if (convergedInternal==false)
				{
					rLogger << "state with fixed yield conditions did not converge. norm Residual " << residual.squaredNorm() << "\n";
					for (int count=0; count<2; count++)
					{
						if (yieldConditionFlag(count)==ACTIVE)
						{
							rLogger << "     yield condition "<< count+1 << " " << yieldCondition(count) << " ("<<toleranceYieldSurface*f_ct << ")"<<"\n";
							rLogger << "     deltaGamma "<< count+1 << " " << deltaGamma(count) << "\n";
						}
						else
							rLogger << "     yield condition "<< count+1 << " " << yieldCondition(count) << "\n";

					}
				}
#endif
			}
		}
		catch (NuTo::MechanicsException& e)
		{
			e.AddMessage("[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] Error performing return mapping procedure.");
			throw e;
		}
		catch (...)
		{
			   throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] Error performing return mapping procedure.");
		}

		if (convergedInternal)
		{
#ifdef ENABLE_DEBUG
			rLogger << "numberOfInternalIterations " << numberOfInternalIterations -  prevNumberOfInternalIterations<< "(" << numberOfInternalIterations << ")" << "\n";
			rLogger << "convergence for external cutback factor" << "\n";
			prevNumberOfInternalIterations = numberOfInternalIterations;
#endif

			//update equivalente plastic strain
			deltaPlasticStrain = rEpsilonP - lastPlastStrain;
			rDeltaEqPlasticStrain += sqrt(deltaPlasticStrain[0]*deltaPlasticStrain(0)+deltaPlasticStrain(1)*deltaPlasticStrain(1)+
							0.5 *(deltaPlasticStrain(2)*deltaPlasticStrain(2))+deltaPlasticStrain(3)*deltaPlasticStrain(3));

//            double tmp(sqrt(rEpsilonP(0)*rEpsilonP(0) + rEpsilonP(1)*rEpsilonP(1)+ 0.5*rEpsilonP(2)*rEpsilonP(2)+ rEpsilonP(3)*rEpsilonP(3)));
//            rLogger << "\n" << "rDeltaEqPlasticStrain " << rDeltaEqPlasticStrain << "Norm of epsilonp " << tmp << "delta between " << rDeltaEqPlasticStrain - tmp << "\n"<< "\n";
//            rLogger << "\n" << "rEpsilonP " << rEpsilonP << "Norm of epsilonp " << "\n";

			if (cutbackFactorExternal==1.)
			{
				convergedExternal=true;
			}
			else
			{
				lastPlastStrain = rEpsilonP;
				lastDeltaEqPlasticStrain = rDeltaEqPlasticStrain;

				if (numberOfInternalIterations<10)
					deltaCutbackFactorExternal*=1.5;
				cutbackFactorExternal +=deltaCutbackFactorExternal;
				if (cutbackFactorExternal>1)
				{
					deltaCutbackFactorExternal -= cutbackFactorExternal -1.;
					cutbackFactorExternal = 1.;
				}
			}
		}
		else
		{
			// no convergence in return mapping
			deltaCutbackFactorExternal*=0.5;
			cutbackFactorExternal-=deltaCutbackFactorExternal;
			//#ifdef ENABLE_DEBUG
						rLogger << "decrease external cutback factor to " << deltaCutbackFactorExternal << "\n";
			//#endif
		}

	}
	if (cutbackFactorExternal<=minCutbackFactor)
	{
		rLogger << "[NuTo::GradientDamagePlasticityEngineeringStress::ReturnMapping2D] No convergence can be obtained in the return mapping procedure." << "n";
		return Error::NO_CONVERGENCE;
	}

    return Error::SUCCESSFUL;
}
//! @brief calculates the first and second derivative of the second Rankine yield surface with respect to the stress
//! @param rStress current stress
//! @param rSigma_1 first principal stress
//! @param rSigma_2 second principal stress
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
//! @return yield condition
double NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine1DRounded(Eigen::Matrix<double,1,1>& rStress, double rFct)const
{
/*    double value_sum  = 0.5*(rStress(0)+rStress(1));
    double help_scalar = (rStress(0)-rStress(1));
    double value_sqrt = sqrt(help_scalar*help_scalar+4.*rStress(2)*rStress(2));
    double sigma_1(value_sum+0.5*value_sqrt);
    double sigma_2(value_sum-0.5*value_sqrt);
    // (rounded) Rankine
#ifdef ENABLE_DEBUG
    std::cout << "p1 " << sigma_1 << " p2 " << sigma_2 << " p3 " << rStress(3) << "\n" << "\n";
#endif
    if (rStress(3)<0)
    {
        //sigma_3 is negative
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            if (rStress(3)>sigma_1)
            {
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " all negative f1" << "\n";
#endif
                return rStress(3)-rFct;
            }
            else
            {
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " all negative f3" << "\n";
#endif
                return sigma_1-rFct;
            }
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                //sigma_2 is negative
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1" << "\n";
#endif
                return sigma_1 - rFct;
            }
            else
            {
                //sigma_2 is positive
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1,f2" << "\n";
#endif
                return sqrt(sigma_1*sigma_1 + sigma_2*sigma_2)-rFct;
            }
        }
    }
    else
    {
        //sigma_3 is positive
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
#ifdef ENABLE_DEBUG
            std::cout << "\n" << " f3" << "\n";
#endif
            return rStress(3)-rFct;
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                //sigma_2 is negative
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1,f3" << "\n";
#endif
                return sqrt(sigma_1*sigma_1 + rStress(3) * rStress(3)) - rFct;
            }
            else
            {
                //sigma_2 is positive
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1,f2,f3" << "\n";
#endif
                return sqrt(sigma_1*sigma_1 + sigma_2*sigma_2 + rStress(3) * rStress(3))-rFct;
            }
        }
    }
    */
	throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine1DRounded] not implemented.");
	return 1000;
}

//! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rStress current stress
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
void NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine1DRoundedDerivatives(Eigen::Matrix<double,1,1>& rdF_dSigma,Eigen::Matrix<double,1,1>* rd2F_d2Sigma,
        Eigen::Matrix<double,1,1>& rStress)const
{
/*    double value_sum  = 0.5*(rStress(0)+rStress(1));
    double help_diff = (rStress(0)-rStress(1));
    double value_sqrt = sqrt(help_diff*help_diff+4.*rStress(2)*rStress(2));
    double sigma_1(value_sum+0.5*value_sqrt);
    double sigma_2(value_sum-0.5*value_sqrt);

    double factor;
    if (rStress(3)<0)
    {
        //sigma_3 is negative
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            if (rStress(3)>sigma_1)
            {
                rdF_dSigma(0) = 0.;
                rdF_dSigma(1) = 0.;
                rdF_dSigma(2) = 0.;
                rdF_dSigma(3) = 1.;
                //second derivatives are not required
            }
            else
            {
                // f = f(sigma_1)
                if (value_sqrt<1e-12)
                    throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

                double yield = sqrt(sigma_1*sigma_1 + rStress(3)*rStress(3));
                factor = 1./yield;
                double dsigma1_dx = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
                double dsigma1_dy = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
                double dsigma1_dxy = 2.*rStress(2)/value_sqrt;

                rdF_dSigma(0) = factor*sigma_1*dsigma1_dx;
                rdF_dSigma(1) = factor*sigma_1*dsigma1_dy;
                rdF_dSigma(2) = factor*sigma_1*dsigma1_dxy;
                rdF_dSigma(3) = factor*rStress(3);
                //second derivatives are not required
            }
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                // f = f(sigma_1)
                if (value_sqrt<1e-12)
                    throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

                rdF_dSigma(0) = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
                rdF_dSigma(1) = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
                rdF_dSigma(2) = 2*rStress(2)/value_sqrt;
                rdF_dSigma(3) = 0.;

                // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                if (rd2F_d2Sigma!=0)
                {
                    double factor = 1./(value_sqrt*value_sqrt*value_sqrt);
                    (*rd2F_d2Sigma)(0,0) = factor*2.*rStress(2)*rStress(2);
                    (*rd2F_d2Sigma)(0,1) = factor*(-2.)*rStress(2)*rStress(2);
                    (*rd2F_d2Sigma)(0,2) = factor*(2.)*(rStress(1)-rStress(0))*rStress(2);
                    (*rd2F_d2Sigma)(0,3) = 0;
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*2.*rStress(2)*rStress(2);
                    (*rd2F_d2Sigma)(1,2) = factor*(2.)*(rStress(0)-rStress(1))*rStress(2);
                    (*rd2F_d2Sigma)(1,3) = 0;
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = factor*(2.)*(rStress(0)-rStress(1))*(rStress(0)-rStress(1));
                    (*rd2F_d2Sigma)(2,3) = 0.;
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = 0.;
                }
            }
            else
            {
                // f = f(sigma_1, sigma_2)
                //sigma_2 is positive
                factor = 1./sqrt(rStress(0)*rStress(0)+rStress(1)*rStress(1)+
                                 2*rStress(2)*rStress(2));
                rdF_dSigma(0) = rStress(0)*factor;
                rdF_dSigma(1) = rStress(1)*factor;
                rdF_dSigma(2) = 2.*rStress(2)*factor;
                rdF_dSigma(3) = 0.;

                if (rd2F_d2Sigma!=0)
                {
                    // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                    factor = sqrt(rStress(0)*rStress(0)+rStress(1)*rStress(1)+2*rStress(2)*rStress(2));
                    factor = 1./(factor*factor*factor);

                    (*rd2F_d2Sigma)(0,0) = factor*(rStress(1)*rStress(1)+2.*rStress(2)*rStress(2));
                    (*rd2F_d2Sigma)(0,1) = -factor*rStress(0)*rStress(1);
                    (*rd2F_d2Sigma)(0,2) = (-2.)*factor*(rStress(0)*rStress(2));
                    (*rd2F_d2Sigma)(0,3) = 0;
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*(rStress(0)*rStress(0)+2.*rStress(2)*rStress(2));
                    (*rd2F_d2Sigma)(1,2) = (-2.)*factor*(rStress(1)*rStress(2));
                    (*rd2F_d2Sigma)(1,3) = 0;
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = (2.)*factor*(rStress(0)*rStress(0)+rStress(1)*rStress(1));
                    (*rd2F_d2Sigma)(2,3) = 0.;
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = 0.;
                }
            }
        }
    }
    else
    {
        //sigma_3 is positive
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            // f = f(sigma_3)
            rdF_dSigma(0) = 0.;
            rdF_dSigma(1) = 0.;
            rdF_dSigma(2) = 0.;
            rdF_dSigma(3) = 1.;

            if (rd2F_d2Sigma!=0)
            {
                (*rd2F_d2Sigma).setZero(4,4);
            }
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                //sigma_2 is negative
                //f = f( sigma_1,sigma_3)

                double yield = sqrt(sigma_1*sigma_1 + rStress(3)*rStress(3));
                factor = 1./yield;
                double dsigma1_dx = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
                double dsigma1_dy = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
                double dsigma1_dxy = 2.*rStress(2)/value_sqrt;

                rdF_dSigma(0) = factor*sigma_1*dsigma1_dx;
                rdF_dSigma(1) = factor*sigma_1*dsigma1_dy;
                rdF_dSigma(2) = factor*sigma_1*dsigma1_dxy;
                rdF_dSigma(3) = factor*rStress(3);

                if (rd2F_d2Sigma!=0)
                {
                    double factor2 = 1./(value_sqrt*value_sqrt*value_sqrt);
                    double dsigma1_dx2    = 2.*rStress(2)*rStress(2)*factor2;
                    double dsigma1_dxdy   = -dsigma1_dx2;
                    double dsigma1_dxdxy  = 2.*rStress(2)*(rStress(1)-rStress(0))*factor2;
                    double dsigma1_dy2    = dsigma1_dx2;
                    double dsigma1_dydxy  = 2.*rStress(2)*(rStress(0)-rStress(1))*factor2;
                    double dsigma1_dxy2   = 2.*help_diff*help_diff*factor2;

//                    rLogger << "d2_sigma_1 d2sigma_i " << dsigma1_dx2 << " " << dsigma1_dxdy << " " << dsigma1_dxdxy << " " << dsigma1_dy2 << " " << dsigma1_dydxy << " " << dsigma1_dxy2 << "\n";

                    factor *= factor;
                    // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                    (*rd2F_d2Sigma)(0,0) = factor*((dsigma1_dx*dsigma1_dx +sigma_1*dsigma1_dx2  )*yield-sigma_1*dsigma1_dx*rdF_dSigma(0));
                    (*rd2F_d2Sigma)(0,1) = factor*((dsigma1_dx*dsigma1_dy +sigma_1*dsigma1_dxdy )*yield-sigma_1*dsigma1_dx*rdF_dSigma(1));
                    (*rd2F_d2Sigma)(0,2) = factor*((dsigma1_dx*dsigma1_dxy+sigma_1*dsigma1_dxdxy)*yield-sigma_1*dsigma1_dx*rdF_dSigma(2));
                    (*rd2F_d2Sigma)(0,3) = -factor*rStress(3)*rdF_dSigma(0);
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*((dsigma1_dy*dsigma1_dy +sigma_1*dsigma1_dy2  )*yield-sigma_1*dsigma1_dy*rdF_dSigma(1));
                    (*rd2F_d2Sigma)(1,2) = factor*((dsigma1_dy*dsigma1_dxy+sigma_1*dsigma1_dydxy)*yield-sigma_1*dsigma1_dy*rdF_dSigma(2));
                    (*rd2F_d2Sigma)(1,3) = -factor*rStress(3)*rdF_dSigma(1);
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = factor*((dsigma1_dxy*dsigma1_dxy +sigma_1*dsigma1_dxy2  )*yield-sigma_1*dsigma1_dxy*rdF_dSigma(2));
                    (*rd2F_d2Sigma)(2,3) = -factor*rStress(3)*rdF_dSigma(2);
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = factor*(yield-rStress(3)*rdF_dSigma(3));

//                    rLogger << "rd2F_d2Sigma " << "\n" << *rd2F_d2Sigma << dsigma1_dxy2 << "\n";
                }
            }
            else
            {
                //sigma_2 is positive
                //f = f( sigma_1,sigma_2, sigma_3)
                factor = 1./sqrt(rStress(0)*rStress(0) + rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
                rdF_dSigma(0) = factor*rStress(0);
                rdF_dSigma(1) = factor*rStress(1);
                rdF_dSigma(2) = factor*2.*rStress(2);
                rdF_dSigma(3) = factor*rStress(3);

                if (rd2F_d2Sigma!=0)
                {

                    factor*=factor*factor;
                    (*rd2F_d2Sigma)(0,0) = factor * (rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
                    (*rd2F_d2Sigma)(0,1) = -factor * rStress(1) * rStress(0);
                    (*rd2F_d2Sigma)(0,2) = -2.*factor * rStress(2) * rStress(0);
                    (*rd2F_d2Sigma)(0,3) = -factor * rStress(3) * rStress(0);
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor * (rStress(0)*rStress(0) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
                    (*rd2F_d2Sigma)(1,2) = -2.*factor * rStress(2) * rStress(1);
                    (*rd2F_d2Sigma)(1,3) = -factor * rStress(3) * rStress(1);
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = 2.*factor*(rStress(0)*rStress(0) + rStress(1)*rStress(1) + rStress(3)*rStress(3));
                    (*rd2F_d2Sigma)(2,3) = -2.*factor * rStress(3) * rStress(2);
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = factor * (rStress(0)*rStress(0) + rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2));
                }
            }
        }
    }

#ifdef ENABLE_DEBUG
    //check the calculation
    double delta(1e-8);
    if (rd2F_d2Sigma!=0 && rStress(0)>1.96)
    {
        //rLogger << "\n" << "check yield surface and derivatives" << "\n";
        //rLogger << "sigmas " << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
        Eigen::Matrix<double,4,1> rdF_dSigma1, rdF_dSigma2,rdF_dSigmaCDF;
        Eigen::Matrix<double,4,4> rd2F_d2SigmaCDF;
        double fct(0.);
        double f1 = YieldSurfaceRankine2DRounded(rStress, fct);
        YieldSurfaceRankine2DRoundedDerivatives(rdF_dSigma1, 0, rStress);
        for (int count=0; count<4; count++)
        {
            rStress(count)+= delta;
            double f2 = YieldSurfaceRankine2DRounded(rStress, fct);
            rdF_dSigmaCDF(count) = (f2-f1)/delta;

            YieldSurfaceRankine2DRoundedDerivatives(rdF_dSigma2, 0, rStress);
            rd2F_d2SigmaCDF.row(count) = (rdF_dSigma2-rdF_dSigma1)/delta;

            rStress(count)-=delta;
        }

        if ((rdF_dSigmaCDF-rdF_dSigma).array().abs().maxCoeff()>1e-1)
        {
            std::cout << "sigmas principal" << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
            std::cout << "sigmas " << rStress(0) << " " << rStress(1) <<  " " << rStress(2) <<  " " <<rStress(3) << "\n";
            std::cout << "error first derivative " << (rdF_dSigmaCDF-rdF_dSigma).array().abs().maxCoeff() << "\n";

            std::cout<< "rdF_dSigma " << "\n" << rdF_dSigma << "\n"<< "\n";
            std::cout<< "rdF_dSigmaCDF " << "\n" << rdF_dSigmaCDF << "\n"<< "\n";
            throw MechanicsException("[NuTo::Mechanics::GradientDamagePlasticityEngineeringStress] Error calculating first derivative of yield function.");
        }

        // fabs is checked, since of the type of the yield surface changes, the second derivatives are likely to change as well
        if (fabs(sigma_1)>delta && fabs(sigma_2)>delta && fabs(rStress(3))>delta)
        {
            if ((rd2F_d2SigmaCDF-(*rd2F_d2Sigma)).array().abs().maxCoeff()>1e-1 && fabs(sigma_1)>delta && fabs(sigma_2)>delta && fabs(rStress(3))>delta)
            {
                std::cout << "sigmas principal" << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
                std::cout << "sigmas " << rStress(0) << " " << rStress(1) <<  " " << rStress(2) <<  " " <<rStress(3) << "\n";
                std::cout << "error second derivatives " << (rd2F_d2SigmaCDF-(*rd2F_d2Sigma)).array().abs().maxCoeff() << "\n";

                std::cout<< "rd2F_d2SigmaCDF " << "\n" << rd2F_d2SigmaCDF << "\n"<< "\n";
                std::cout<< "rd2F_d2Sigma " << "\n" << (*rd2F_d2Sigma) << "\n"<< "\n";
                throw MechanicsException("[NuTo::Mechanics::GradientDamagePlasticityEngineeringStress] Error calculating second derivative of yield function.");
            }
        }
        else
        {
            std::cout << "at least one principal stress close to zero, stiffness matrix is not reliable" << "\n";
        }
    }
#endif
*/
}



//! @brief calculates the first and second derivative of the second Rankine yield surface with respect to the stress
//! @param rStress current stress
//! @param rSigma_1 first principal stress
//! @param rSigma_2 second principal stress
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
//! @return yield condition
double NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded(Eigen::Matrix<double,4,1>& rStress, double rFct)const
{
    double value_sum  = 0.5*(rStress(0)+rStress(1));
    double help_scalar = (rStress(0)-rStress(1));
    double value_sqrt = sqrt(help_scalar*help_scalar+4.*rStress(2)*rStress(2));
    double sigma_1(value_sum+0.5*value_sqrt);
    double sigma_2(value_sum-0.5*value_sqrt);
    //* (rounded) Rankine
#ifdef ENABLE_DEBUG
    std::cout << "p1 " << sigma_1 << " p2 " << sigma_2 << " p3 " << rStress(3) << "\n" << "\n";
#endif
    if (rStress(3)<0)
    {
        //sigma_3 is negative
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            if (rStress(3)>sigma_1)
            {
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " all negative f1" << "\n";
#endif
                return rStress(3)-rFct;
            }
            else
            {
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " all negative f3" << "\n";
#endif
                return sigma_1-rFct;
            }
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                //sigma_2 is negative
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1" << "\n";
#endif
                return sigma_1 - rFct;
            }
            else
            {
                //sigma_2 is positive
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1,f2" << "\n";
#endif
                return sqrt(sigma_1*sigma_1 + sigma_2*sigma_2)-rFct;
            }
        }
    }
    else
    {
        //sigma_3 is positive
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
#ifdef ENABLE_DEBUG
            std::cout << "\n" << " f3" << "\n";
#endif
            return rStress(3)-rFct;
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                //sigma_2 is negative
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1,f3" << "\n";
#endif
                return sqrt(sigma_1*sigma_1 + rStress(3) * rStress(3)) - rFct;
            }
            else
            {
                //sigma_2 is positive
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1,f2,f3" << "\n";
#endif
                return sqrt(sigma_1*sigma_1 + sigma_2*sigma_2 + rStress(3) * rStress(3))-rFct;
            }
        }
    }
}

//! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rStress current stress
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
void NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRoundedDerivatives(Eigen::Matrix<double,4,1>& rdF_dSigma,Eigen::Matrix<double,4,4>* rd2F_d2Sigma,
        Eigen::Matrix<double,4,1>& rStress)const
{
    double value_sum  = 0.5*(rStress(0)+rStress(1));
    double help_diff = (rStress(0)-rStress(1));
    double value_sqrt = sqrt(help_diff*help_diff+4.*rStress(2)*rStress(2));
    double sigma_1(value_sum+0.5*value_sqrt);
    double sigma_2(value_sum-0.5*value_sqrt);

    double factor;
    if (rStress(3)<0)
    {
        //sigma_3 is negative
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            if (rStress(3)>sigma_1)
            {
                rdF_dSigma(0) = 0.;
                rdF_dSigma(1) = 0.;
                rdF_dSigma(2) = 0.;
                rdF_dSigma(3) = 1.;
                //second derivatives are not required
            }
            else
            {
                // f = f(sigma_1)
                if (value_sqrt<1e-12)
                    throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

                double yield = sqrt(sigma_1*sigma_1 + rStress(3)*rStress(3));
                factor = 1./yield;
                double dsigma1_dx = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
                double dsigma1_dy = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
                double dsigma1_dxy = 2.*rStress(2)/value_sqrt;

                rdF_dSigma(0) = factor*sigma_1*dsigma1_dx;
                rdF_dSigma(1) = factor*sigma_1*dsigma1_dy;
                rdF_dSigma(2) = factor*sigma_1*dsigma1_dxy;
                rdF_dSigma(3) = factor*rStress(3);
                //second derivatives are not required
            }
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                // f = f(sigma_1)
                if (value_sqrt<1e-12)
                    throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

                rdF_dSigma(0) = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
                rdF_dSigma(1) = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
                rdF_dSigma(2) = 2*rStress(2)/value_sqrt;
                rdF_dSigma(3) = 0.;

                // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                if (rd2F_d2Sigma!=0)
                {
                    double factor = 1./(value_sqrt*value_sqrt*value_sqrt);
                    (*rd2F_d2Sigma)(0,0) = factor*2.*rStress(2)*rStress(2);
                    (*rd2F_d2Sigma)(0,1) = factor*(-2.)*rStress(2)*rStress(2);
                    (*rd2F_d2Sigma)(0,2) = factor*(2.)*(rStress(1)-rStress(0))*rStress(2);
                    (*rd2F_d2Sigma)(0,3) = 0;
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*2.*rStress(2)*rStress(2);
                    (*rd2F_d2Sigma)(1,2) = factor*(2.)*(rStress(0)-rStress(1))*rStress(2);
                    (*rd2F_d2Sigma)(1,3) = 0;
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = factor*(2.)*(rStress(0)-rStress(1))*(rStress(0)-rStress(1));
                    (*rd2F_d2Sigma)(2,3) = 0.;
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = 0.;
                }
            }
            else
            {
                // f = f(sigma_1, sigma_2)
                //sigma_2 is positive
                factor = 1./sqrt(rStress(0)*rStress(0)+rStress(1)*rStress(1)+
                                 2*rStress(2)*rStress(2));
                rdF_dSigma(0) = rStress(0)*factor;
                rdF_dSigma(1) = rStress(1)*factor;
                rdF_dSigma(2) = 2.*rStress(2)*factor;
                rdF_dSigma(3) = 0.;

                if (rd2F_d2Sigma!=0)
                {
                    // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                    factor = sqrt(rStress(0)*rStress(0)+rStress(1)*rStress(1)+2*rStress(2)*rStress(2));
                    factor = 1./(factor*factor*factor);

                    (*rd2F_d2Sigma)(0,0) = factor*(rStress(1)*rStress(1)+2.*rStress(2)*rStress(2));
                    (*rd2F_d2Sigma)(0,1) = -factor*rStress(0)*rStress(1);
                    (*rd2F_d2Sigma)(0,2) = (-2.)*factor*(rStress(0)*rStress(2));
                    (*rd2F_d2Sigma)(0,3) = 0;
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*(rStress(0)*rStress(0)+2.*rStress(2)*rStress(2));
                    (*rd2F_d2Sigma)(1,2) = (-2.)*factor*(rStress(1)*rStress(2));
                    (*rd2F_d2Sigma)(1,3) = 0;
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = (2.)*factor*(rStress(0)*rStress(0)+rStress(1)*rStress(1));
                    (*rd2F_d2Sigma)(2,3) = 0.;
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = 0.;
                }
            }
        }
    }
    else
    {
        //sigma_3 is positive
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            // f = f(sigma_3)
            rdF_dSigma(0) = 0.;
            rdF_dSigma(1) = 0.;
            rdF_dSigma(2) = 0.;
            rdF_dSigma(3) = 1.;

            if (rd2F_d2Sigma!=0)
            {
                (*rd2F_d2Sigma).setZero(4,4);
            }
        }
        else
        {
            //sigma_1 is positive
            if (sigma_2<0)
            {
                //sigma_2 is negative
                //f = f( sigma_1,sigma_3)

                double yield = sqrt(sigma_1*sigma_1 + rStress(3)*rStress(3));
                factor = 1./yield;
                double dsigma1_dx = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
                double dsigma1_dy = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
                double dsigma1_dxy = 2.*rStress(2)/value_sqrt;

                rdF_dSigma(0) = factor*sigma_1*dsigma1_dx;
                rdF_dSigma(1) = factor*sigma_1*dsigma1_dy;
                rdF_dSigma(2) = factor*sigma_1*dsigma1_dxy;
                rdF_dSigma(3) = factor*rStress(3);

                if (rd2F_d2Sigma!=0)
                {
/*
                    rLogger << "active yieldsurfaces : p1, p3 " << rd2F_d2Sigma << "\n";
                    rLogger << "yield " << yield << "\n";
                    rLogger << "d_sigma_1 dsigma_i " << dsigma1_dx << " " << dsigma1_dy << " " << dsigma1_dxy << "\n";
                    rLogger << "rdF_dSigma " << "\n" << rdF_dSigma << "\n";
*/
                    double factor2 = 1./(value_sqrt*value_sqrt*value_sqrt);
                    double dsigma1_dx2    = 2.*rStress(2)*rStress(2)*factor2;
                    double dsigma1_dxdy   = -dsigma1_dx2;
                    double dsigma1_dxdxy  = 2.*rStress(2)*(rStress(1)-rStress(0))*factor2;
                    double dsigma1_dy2    = dsigma1_dx2;
                    double dsigma1_dydxy  = 2.*rStress(2)*(rStress(0)-rStress(1))*factor2;
                    double dsigma1_dxy2   = 2.*help_diff*help_diff*factor2;

//                    rLogger << "d2_sigma_1 d2sigma_i " << dsigma1_dx2 << " " << dsigma1_dxdy << " " << dsigma1_dxdxy << " " << dsigma1_dy2 << " " << dsigma1_dydxy << " " << dsigma1_dxy2 << "\n";

                    factor *= factor;
                    // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                    (*rd2F_d2Sigma)(0,0) = factor*((dsigma1_dx*dsigma1_dx +sigma_1*dsigma1_dx2  )*yield-sigma_1*dsigma1_dx*rdF_dSigma(0));
                    (*rd2F_d2Sigma)(0,1) = factor*((dsigma1_dx*dsigma1_dy +sigma_1*dsigma1_dxdy )*yield-sigma_1*dsigma1_dx*rdF_dSigma(1));
                    (*rd2F_d2Sigma)(0,2) = factor*((dsigma1_dx*dsigma1_dxy+sigma_1*dsigma1_dxdxy)*yield-sigma_1*dsigma1_dx*rdF_dSigma(2));
                    (*rd2F_d2Sigma)(0,3) = -factor*rStress(3)*rdF_dSigma(0);
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*((dsigma1_dy*dsigma1_dy +sigma_1*dsigma1_dy2  )*yield-sigma_1*dsigma1_dy*rdF_dSigma(1));
                    (*rd2F_d2Sigma)(1,2) = factor*((dsigma1_dy*dsigma1_dxy+sigma_1*dsigma1_dydxy)*yield-sigma_1*dsigma1_dy*rdF_dSigma(2));
                    (*rd2F_d2Sigma)(1,3) = -factor*rStress(3)*rdF_dSigma(1);
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = factor*((dsigma1_dxy*dsigma1_dxy +sigma_1*dsigma1_dxy2  )*yield-sigma_1*dsigma1_dxy*rdF_dSigma(2));
                    (*rd2F_d2Sigma)(2,3) = -factor*rStress(3)*rdF_dSigma(2);
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = factor*(yield-rStress(3)*rdF_dSigma(3));

//                    rLogger << "rd2F_d2Sigma " << "\n" << *rd2F_d2Sigma << dsigma1_dxy2 << "\n";
                }
            }
            else
            {
                //sigma_2 is positive
                //f = f( sigma_1,sigma_2, sigma_3)
                factor = 1./sqrt(rStress(0)*rStress(0) + rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
                rdF_dSigma(0) = factor*rStress(0);
                rdF_dSigma(1) = factor*rStress(1);
                rdF_dSigma(2) = factor*2.*rStress(2);
                rdF_dSigma(3) = factor*rStress(3);

                if (rd2F_d2Sigma!=0)
                {

                    factor*=factor*factor;
                    (*rd2F_d2Sigma)(0,0) = factor * (rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
                    (*rd2F_d2Sigma)(0,1) = -factor * rStress(1) * rStress(0);
                    (*rd2F_d2Sigma)(0,2) = -2.*factor * rStress(2) * rStress(0);
                    (*rd2F_d2Sigma)(0,3) = -factor * rStress(3) * rStress(0);
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor * (rStress(0)*rStress(0) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
                    (*rd2F_d2Sigma)(1,2) = -2.*factor * rStress(2) * rStress(1);
                    (*rd2F_d2Sigma)(1,3) = -factor * rStress(3) * rStress(1);
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = 2.*factor*(rStress(0)*rStress(0) + rStress(1)*rStress(1) + rStress(3)*rStress(3));
                    (*rd2F_d2Sigma)(2,3) = -2.*factor * rStress(3) * rStress(2);
                    (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
                    (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
                    (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
                    (*rd2F_d2Sigma)(3,3) = factor * (rStress(0)*rStress(0) + rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2));
                }
            }
        }
    }

#ifdef ENABLE_DEBUG
    //check the calculation
    double delta(1e-8);
    if (rd2F_d2Sigma!=0 && rStress(0)>1.96)
    {
        //rLogger << "\n" << "check yield surface and derivatives" << "\n";
        //rLogger << "sigmas " << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
        Eigen::Matrix<double,4,1> rdF_dSigma1, rdF_dSigma2,rdF_dSigmaCDF;
        Eigen::Matrix<double,4,4> rd2F_d2SigmaCDF;
        double fct(0.);
        double f1 = YieldSurfaceRankine2DRounded(rStress, fct);
        YieldSurfaceRankine2DRoundedDerivatives(rdF_dSigma1, 0, rStress);
        for (int count=0; count<4; count++)
        {
            rStress(count)+= delta;
            double f2 = YieldSurfaceRankine2DRounded(rStress, fct);
            rdF_dSigmaCDF(count) = (f2-f1)/delta;

            YieldSurfaceRankine2DRoundedDerivatives(rdF_dSigma2, 0, rStress);
            rd2F_d2SigmaCDF.row(count) = (rdF_dSigma2-rdF_dSigma1)/delta;

            rStress(count)-=delta;
        }

        if ((rdF_dSigmaCDF-rdF_dSigma).array().abs().maxCoeff()>1e-1)
        {
            std::cout << "sigmas principal" << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
            std::cout << "sigmas " << rStress(0) << " " << rStress(1) <<  " " << rStress(2) <<  " " <<rStress(3) << "\n";
            std::cout << "error first derivative " << (rdF_dSigmaCDF-rdF_dSigma).array().abs().maxCoeff() << "\n";

            std::cout<< "rdF_dSigma " << "\n" << rdF_dSigma << "\n"<< "\n";
            std::cout<< "rdF_dSigmaCDF " << "\n" << rdF_dSigmaCDF << "\n"<< "\n";
            throw MechanicsException("[NuTo::Mechanics::GradientDamagePlasticityEngineeringStress] Error calculating first derivative of yield function.");
        }

        // fabs is checked, since of the type of the yield surface changes, the second derivatives are likely to change as well
        if (fabs(sigma_1)>delta && fabs(sigma_2)>delta && fabs(rStress(3))>delta)
        {
            if ((rd2F_d2SigmaCDF-(*rd2F_d2Sigma)).array().abs().maxCoeff()>1e-1 && fabs(sigma_1)>delta && fabs(sigma_2)>delta && fabs(rStress(3))>delta)
            {
                std::cout << "sigmas principal" << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
                std::cout << "sigmas " << rStress(0) << " " << rStress(1) <<  " " << rStress(2) <<  " " <<rStress(3) << "\n";
                std::cout << "error second derivatives " << (rd2F_d2SigmaCDF-(*rd2F_d2Sigma)).array().abs().maxCoeff() << "\n";

                std::cout<< "rd2F_d2SigmaCDF " << "\n" << rd2F_d2SigmaCDF << "\n"<< "\n";
                std::cout<< "rd2F_d2Sigma " << "\n" << (*rd2F_d2Sigma) << "\n"<< "\n";
                throw MechanicsException("[NuTo::Mechanics::GradientDamagePlasticityEngineeringStress] Error calculating second derivative of yield function.");
            }
        }
        else
        {
            std::cout << "at least one principal stress close to zero, stiffness matrix is not reliable" << "\n";
        }
    }
#endif
}

//! @brief calculates the first and second derivative of the second Rankine yield surface with respect to the stress
//! @param rStress current stress
//! @param rBeta parameter of the Drucker-Prager yield surface
//! @param rHP parameter of the Drucker-Prager yield surface
//! @return yield condition
double NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceDruckerPrager2D(Eigen::Matrix<double,4,1>& rStress, double rBeta, double rHP)const
{
    double invariante_1 = rStress(0)+rStress(1)+rStress(3);

    //*******************************************************************
    //*    second invariante                                                 *
    //*******************************************************************
    double invariante_2 = ((rStress(0)-rStress(1))*(rStress(0)-rStress(1))+
                            (rStress(1)-rStress(3))*(rStress(1)-rStress(3))+
                            (rStress(0)-rStress(3))*(rStress(0)-rStress(3)))/6.+
                             rStress(2)*rStress(2);

    //*******************************************************************
    //*    F_DEV:    term of stress-deviator with/without kinematic hard.    *
    //*******************************************************************
    double F_DEV = sqrt( invariante_2 ) ;

    //*******************************************************************
    //*    check yield first time                                            *
    //*******************************************************************
    double F_BETA = invariante_1 * rBeta/3. ;
    double F_FLOW = -rHP;

    //* calculate yield condition and yield condition flag (active or not)
    //* Drucker Prager
    return F_BETA + F_DEV + F_FLOW;
}

//! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rStress current stress
//! @param rBETA parameter of the Drucker Prager yield surface
//! @return false, if the stress is on the hydrostatic axis otherwise true
bool NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceDruckerPrager2DDerivatives(Eigen::Matrix<double,4,1>& rdF_dSigma,Eigen::Matrix<double,4,4>* rd2F_d2Sigma,
        Eigen::Matrix<double,4,1>& rStress, double rBETA)const
{
    //*******************************************************************
    //*    second invariante                                                 *
    //*******************************************************************
    double invariante_2 = ((rStress(0)-rStress(1))*(rStress(0)-rStress(1))+
                            (rStress(1)-rStress(3))*(rStress(1)-rStress(3))+
                            (rStress(0)-rStress(3))*(rStress(0)-rStress(3)))/6.+
                             rStress(2)*rStress(2);

    if (fabs(invariante_2)<1e-12)
        return false;
    double factor = 1./(sqrt(invariante_2)*6.);

    rdF_dSigma(0) = factor * (2.*rStress(0)-rStress(1)-rStress(3))+rBETA/3.;
    rdF_dSigma(1) = factor * (2.*rStress(1)-rStress(0)-rStress(3))+rBETA/3.;
    rdF_dSigma(2) = factor * (6.*rStress(2)); /* vector notation from second order tensor */
    rdF_dSigma(3) = factor * (2.*rStress(3)-rStress(0)-rStress(1))+rBETA/3.;

    if (rd2F_d2Sigma!=0)
    {
        factor = 1./(invariante_2*sqrt(invariante_2)*12.);
        (*rd2F_d2Sigma)(0,0) = factor * (rStress(1)*rStress(1)-2.*rStress(1)*rStress(3)+rStress(3)*rStress(3)+4.*rStress(2)*rStress(2));
        (*rd2F_d2Sigma)(0,1) = factor * (-rStress(0)*rStress(1)+rStress(0)*rStress(3)+rStress(1)*rStress(3)-rStress(3)*rStress(3)-2.*rStress(2)*rStress(2));
        (*rd2F_d2Sigma)(0,2) = factor * (-2.*rStress(2)*(2.*rStress(0)-rStress(1)-rStress(3)));
        (*rd2F_d2Sigma)(0,3) = factor * (rStress(0)*rStress(1)-rStress(0)*rStress(3)-rStress(1)*rStress(1)+rStress(1)*rStress(3)-2.*rStress(2)*rStress(2));
        (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
        (*rd2F_d2Sigma)(1,1) = factor * (rStress(0)*rStress(0)-2.*rStress(0)*rStress(3)+rStress(3)*rStress(3)+4.*rStress(2)*rStress(2));
        (*rd2F_d2Sigma)(1,2) = factor * (2.*rStress(2)*(rStress(0)-2.*rStress(1)+rStress(3)));
        (*rd2F_d2Sigma)(1,3) = factor * (-rStress(0)*rStress(0)+rStress(0)*rStress(1)+rStress(0)*rStress(3)-rStress(1)*rStress(3)-2.*rStress(2)*rStress(2));
        (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
        (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
        (*rd2F_d2Sigma)(2,2) = factor * 4.*(rStress(0)*rStress(0)-rStress(0)*rStress(1)-rStress(0)*rStress(3)+
                                        rStress(1)*rStress(1)-rStress(1)*rStress(3)+rStress(3)*rStress(3));
        (*rd2F_d2Sigma)(2,3) = factor * 2. *(rStress(0)+rStress(1)-2.*rStress(3))*rStress(2);
        (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
        (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
        (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
        (*rd2F_d2Sigma)(3,3) = factor * (rStress(0)*rStress(0)-2*rStress(0)*rStress(1)+rStress(1)*rStress(1)+4.*rStress(2)*rStress(2));
    }
    return true;
}

#define tolerance_D 1e-10
#define tolerance_q sqrt(tolerance_D)
//! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
//! @param rStress current stress
//! @param rFct return tensile strength
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @return yield function
double NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine3DRounded(
		const Eigen::Matrix<double,6,1>& rStress,
		double rFct,
		Eigen::Matrix<double,6,1>* rdF_dSigma,
		Eigen::Matrix<double,6,6>* rd2F_d2Sigma)const
{
    double F,
    a,
    a2,
    a3,
    b,
    c,
    q,
    p,
    D,
    P=1e99,
    beta=1e99,
    principal[3],
    sqrt_minus_p=1e99,
    sqrt_1_minus_q2_div_P6=1e99,
    q2,
    factor1,
    factor2,
    factor3,
    factor4,
    factor5,
    P2,
    P3,
    P6,
    help_scalar,
    cos_beta=1e99,
    cos_beta_add=1e99,
    cos_beta_sub=1e99,
    sin_beta,
    sin_beta_add,
    sin_beta_sub,
    da_ds[6]={0},
    db_ds[6]={0},
    dc_ds[6]={0},
    dp_ds[6]={0},
    dq_ds[6]={0},
    dP_ds[6]={0},
    dbeta_ds[6],
    db_ds2[36]={0},
    dc_ds2[36]={0},
    dp_ds2[36]={0},
    dq_ds2[36]={0},
    dP_ds2[36]={0},
    dbeta_ds2[36]={0};

    int num_pos; //number of positive eigenvalues

    if (rd2F_d2Sigma!=0 && rdF_dSigma==0)
    	throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine3DRounded] second derivative can only be calculated with the first derivative.");
    double *dF_dsigma(0),*d2F_d2Sigma(0);
    if (rdF_dSigma!=0)
    	dF_dsigma =  (*rdF_dSigma).data();
    if (rd2F_d2Sigma!=0)
    	d2F_d2Sigma =  (*rd2F_d2Sigma).data();
    // ***************************
    // 	calculate eigenvalues    *
    // ***************************
    a = -rStress(0)-rStress(1)-rStress(2);
    b =  rStress(0)*rStress(2)
         +rStress(1)*rStress(2)
         +rStress(1)*rStress(0)
         -rStress(3)*rStress(3)
         -rStress(4)*rStress(4)
         -rStress(5)*rStress(5);
    c =  rStress(0)*rStress(4)*rStress(4)
         +rStress(1)*rStress(5)*rStress(5)
         +rStress(2)*rStress(3)*rStress(3)
         -rStress(0)*rStress(1)*rStress(2)
         -2.*rStress(3)*rStress(4)*rStress(5);

    a2 = a*a;
    a3 = a2*a;

    q = a3/27.-a*b/6.+0.5*c;
    p = (b-a2/3.)/3.;
    q2 = q*q;
    D = p*p*p+q2;

    // complex eigenvalues found, but eigenvalues of a symmetric matrix are always real
    //std::cout << "D=" << D << std::endl;
//    std::cout << "stress : " << rStress.transpose() << std::endl ;
    assert(D<tolerance_D);

	if (fabs(q)<tolerance_q)
	{
		//three identical eigenvalues
	    assert(D>-tolerance_D);
		principal[0] = -a/3.;
		principal[1] = principal[0];
		principal[2] = principal[0];
	}
    else
    {
        // either three different real roots or two real roots (one is twice)
    	sqrt_minus_p = sqrt(-p);
        P = q<0 ? -sqrt_minus_p : sqrt_minus_p;
        double qdiffP3 = q/(P*P*P);
        if (qdiffP3>1.)
        {
        	assert(qdiffP3-1<tolerance_D);
        	qdiffP3=1;
        }
        beta = acos(qdiffP3)/3.;
        cos_beta = cos(beta);
        cos_beta_add = cos(beta+M_PI/3.);
        cos_beta_sub = cos(beta-M_PI/3.);
        principal[0] = -2.*P*cos_beta-a/3.;
        principal[1] = 2.*P*cos_beta_add-a/3.;
        principal[2] = 2.*P*cos_beta_sub-a/3.;
    }

    //std::cout << "principal stresses=" <<  principal[0] << " " << principal[1] << " " << principal[2] << std::endl;
	num_pos = (principal[0]>0 ? 1 : 0) + (principal[1]>0 ? 1 : 0) + (principal[2]>0 ? 1 : 0);

	switch (num_pos)
	{
	case 0:
		F = -rFct;
		if (rdF_dSigma!=0)
		{
			(*rdF_dSigma).setZero(6,1);
		}
		if (rd2F_d2Sigma!=0)
		{
			(*rd2F_d2Sigma).setZero(36,1);
		}
		return F;
		break;
	case 1:
	case 2:
		if (rdF_dSigma!=0)
		{
			// a
			da_ds[0] = -1.;
			da_ds[1] = -1.;
			da_ds[2] = -1.;
			da_ds[3] =  0.;
			da_ds[4] =  0.;
			da_ds[5] =  0.;

			// b
			db_ds[0] = rStress(1)+rStress(2);
			db_ds[1] = rStress(0)+rStress(2);
			db_ds[2] = rStress(0)+rStress(1);
			db_ds[3] = -2.*rStress(3);
			db_ds[4] = -2.*rStress(4);
			db_ds[5] = -2.*rStress(5);

			db_ds2[1]  = 1.;
			db_ds2[2]  = 1.;
			db_ds2[8]  = 1.;
			db_ds2[21]  = -2.;
			db_ds2[28] = -2.;
			db_ds2[35] = -2.;

			// c
			dc_ds[0] = rStress(4)*rStress(4)-rStress(1)*rStress(2);
			dc_ds[1] = rStress(5)*rStress(5)-rStress(0)*rStress(2);
			dc_ds[2] = rStress(3)*rStress(3)-rStress(0)*rStress(1);
			dc_ds[3] = 2.*(rStress(2)*rStress(3)-rStress(4)*rStress(5));
			dc_ds[4] = 2.*(rStress(0)*rStress(4)-rStress(3)*rStress(5));
			dc_ds[5] = 2.*(rStress(1)*rStress(5)-rStress(3)*rStress(4));

			dc_ds2[1]  = -rStress(2);
			dc_ds2[2]  = -rStress(1);
			dc_ds2[4]  = 2.*rStress(4);
			dc_ds2[8]  = -rStress(0);
			dc_ds2[11] = 2.*rStress(5);
			dc_ds2[15] = 2.*rStress(3);
			dc_ds2[21] = 2.*rStress(2);
			dc_ds2[22] = -2.*rStress(5);
			dc_ds2[23] = -2.*rStress(4);
			dc_ds2[28] = 2.*rStress(0);
			dc_ds2[29] = -2.*rStress(3);
			dc_ds2[35] = 2.*rStress(1);

			// p
			factor3 = -2./9.;
			factor1 = factor3*a;factor2 = 1./3.;
			dp_ds[0] = factor1*da_ds[0]+factor2*db_ds[0];
			dp_ds[1] = factor1*da_ds[1]+factor2*db_ds[1];
			dp_ds[2] = factor1*da_ds[2]+factor2*db_ds[2];
			dp_ds[3] = factor1*da_ds[3]+factor2*db_ds[3];
			dp_ds[4] = factor1*da_ds[4]+factor2*db_ds[4];
			dp_ds[5] = factor1*da_ds[5]+factor2*db_ds[5];

			if (rd2F_d2Sigma!=0)
			{
				for (int j=0; j<6;j++)
				{
					for (int i=j; i<6;i++)
					{
						dp_ds2[6*j+i] = factor3*da_ds[i]*da_ds[j]+factor2*db_ds2[6*j+i];
					}
				}
			}

			// q
			factor1 = a2/9.-b/6.;factor2 = a/(-6.);
			dq_ds[0] = factor1*da_ds[0]+factor2*db_ds[0]+0.5*dc_ds[0];
			dq_ds[1] = factor1*da_ds[1]+factor2*db_ds[1]+0.5*dc_ds[1];
			dq_ds[2] = factor1*da_ds[2]+factor2*db_ds[2]+0.5*dc_ds[2];
			dq_ds[3] = factor1*da_ds[3]+factor2*db_ds[3]+0.5*dc_ds[3];
			dq_ds[4] = factor1*da_ds[4]+factor2*db_ds[4]+0.5*dc_ds[4];
			dq_ds[5] = factor1*da_ds[5]+factor2*db_ds[5]+0.5*dc_ds[5];

			if (rd2F_d2Sigma!=0)
			{
				factor3 = 2./9.*a;factor4 = -1./6.;
				for (int j=0; j<6;j++)
				{
					for (int i=j; i<6;i++)
					{
						dq_ds2[6*j+i] = (factor3*da_ds[j]+factor4*db_ds[j]) * da_ds[i]
										+factor4*da_ds[j]*db_ds[i]
										+factor2*db_ds2[6*j+i]+
										+0.5*dc_ds2[6*j+i];
					}
				}
			}
			if (fabs(q)>=tolerance_q)
			{
				assert(fabs(p)>tolerance_D);
				// P
				factor1 = q<0 ? 1./(2.*sqrt_minus_p) : -1./(2.*sqrt_minus_p);
				dP_ds[0] = factor1*dp_ds[0];
				dP_ds[1] = factor1*dp_ds[1];
				dP_ds[2] = factor1*dp_ds[2];
				dP_ds[3] = factor1*dp_ds[3];
				dP_ds[4] = factor1*dp_ds[4];
				dP_ds[5] = factor1*dp_ds[5];

				if (rd2F_d2Sigma!=0)
				{
					factor2 = q<0 ? 1./(4.*sqrt_minus_p*(-p)) : -1./(4.*sqrt_minus_p*(-p));
					for (int j=0; j<6;j++)
					{
						for (int i=j; i<6;i++)
						{
							dP_ds2[6*j+i] = factor2*dp_ds[j] * dp_ds[i]+factor1*dp_ds2[6*j+i];
						}
					}
				}
		    }

			if (D<=-tolerance_D)
			{
				// beta
				P2 = P*P;
				P3 = P2*P;
				P6 = P3*P3;
				help_scalar = q/(P3);
				sqrt_1_minus_q2_div_P6 = sqrt(1.-help_scalar*help_scalar);
				factor1 = -1./(3*P3*sqrt_1_minus_q2_div_P6);
				factor2 = q/(P2*P2*sqrt_1_minus_q2_div_P6);
				factor3 = P2/((P6-q2)*sqrt_1_minus_q2_div_P6);
				factor4 = -q*(4.*P6-q2)/((P6-q2)*P2*P3*sqrt_1_minus_q2_div_P6);

				help_scalar = P3*sqrt_1_minus_q2_div_P6;
				factor5 = -q/(3*(help_scalar*help_scalar*help_scalar));

				dbeta_ds[0] = factor1*dq_ds[0]+factor2*dP_ds[0];
				dbeta_ds[1] = factor1*dq_ds[1]+factor2*dP_ds[1];
				dbeta_ds[2] = factor1*dq_ds[2]+factor2*dP_ds[2];
				dbeta_ds[3] = factor1*dq_ds[3]+factor2*dP_ds[3];
				dbeta_ds[4] = factor1*dq_ds[4]+factor2*dP_ds[4];
				dbeta_ds[5] = factor1*dq_ds[5]+factor2*dP_ds[5];

				if (rd2F_d2Sigma!=0)
				{
					for (int j=0; j<6;j++)
					{
						for (int i=j; i<6;i++)
						{
							dbeta_ds2[6*j+i] = (factor5*dq_ds[j]+factor3*dP_ds[j]) * dq_ds[i] + factor1*dq_ds2[6*j+i]
											   +(factor4*dP_ds[j]+factor3*dq_ds[j]) * dP_ds[i] + factor2*dP_ds2[6*j+i];
						}
					}
				}
			} // if (D<=-tolerance_D)
		}// if(calc_derivatives)

		if (num_pos==1)
		{
			if (D<=-tolerance_D)
			{
				//all solutions are real and different
				if (principal[0]>0)
				{

					F = principal[0]-rFct;
					if (rdF_dSigma!=0)
					{
						// calculate derivatives
						sin_beta = sin(beta);
						factor4 = 2.*sin_beta;
						factor1 = P*factor4;
						factor2 = -2.*cos_beta;
						factor3 = -1./3.;
						factor5 = -P*factor2;
						dF_dsigma[0] = factor1*dbeta_ds[0]+factor2*dP_ds[0]+factor3*da_ds[0];
						dF_dsigma[1] = factor1*dbeta_ds[1]+factor2*dP_ds[1]+factor3*da_ds[1];
						dF_dsigma[2] = factor1*dbeta_ds[2]+factor2*dP_ds[2]+factor3*da_ds[2];
						dF_dsigma[3] = factor1*dbeta_ds[3]+factor2*dP_ds[3]+factor3*da_ds[3];
						dF_dsigma[4] = factor1*dbeta_ds[4]+factor2*dP_ds[4]+factor3*da_ds[4];
						dF_dsigma[5] = factor1*dbeta_ds[5]+factor2*dP_ds[5]+factor3*da_ds[5];

						if (rd2F_d2Sigma!=0)
						{
							for (int j=0; j<6;j++)
							{
								for (int i=j; i<6;i++)
								{
									double tmp = (factor5*dbeta_ds[j]+factor4*dP_ds[j]) * dbeta_ds[i] + factor1*dbeta_ds2[6*j+i]
														 +factor4*dbeta_ds[j]*dP_ds[i] + factor2*dP_ds2[6*j+i];
									d2F_d2Sigma[6*j+i] = tmp;
									if (i!=j)
									{
										d2F_d2Sigma[6*i+j] = tmp;
									}
								}
							}
						}
					}
					return F;
				}
				if (principal[1]>0)
				{

					F = principal[1]-rFct;
					if (rdF_dSigma!=0)
					{
						sin_beta_add = sin(beta+M_PI/3.);
						factor4 = -2.*sin_beta_add;
						factor1 = P*factor4;
						factor2 = 2.*cos_beta_add;
						factor3 = -1./3.;
						factor5 = -P*factor2;
						dF_dsigma[0] = factor1*dbeta_ds[0]+factor2*dP_ds[0]+factor3*da_ds[0];
						dF_dsigma[1] = factor1*dbeta_ds[1]+factor2*dP_ds[1]+factor3*da_ds[1];
						dF_dsigma[2] = factor1*dbeta_ds[2]+factor2*dP_ds[2]+factor3*da_ds[2];
						dF_dsigma[3] = factor1*dbeta_ds[3]+factor2*dP_ds[3]+factor3*da_ds[3];
						dF_dsigma[4] = factor1*dbeta_ds[4]+factor2*dP_ds[4]+factor3*da_ds[4];
						dF_dsigma[5] = factor1*dbeta_ds[5]+factor2*dP_ds[5]+factor3*da_ds[5];

						if (rd2F_d2Sigma!=0)
						{
							for (int j=0; j<6;j++)
							{
								for (int i=j; i<6;i++)
								{
									double tmp = (factor5*dbeta_ds[j]+factor4*dP_ds[j]) * dbeta_ds[i] + factor1*dbeta_ds2[6*j+i]
														 +factor4*dbeta_ds[j]*dP_ds[i] + factor2*dP_ds2[6*j+i];
									d2F_d2Sigma[6*j+i] = tmp;
									if (i!=j)
									{
										d2F_d2Sigma[6*i+j] = tmp;
									}
								}
							}
						}
					}
					return F;
				}
				if (principal[2]>0)
				{
					F = principal[2]-rFct;
					if (rdF_dSigma!=0)
					{
						sin_beta_sub = sin(beta-M_PI/3.);
						factor4 = -2.*sin_beta_sub;
						factor1 = P*factor4;
						factor2 = 2.*cos_beta_sub;
						factor3 = -1./3.;
						factor5 = -P*factor2;
						dF_dsigma[0] = factor1*dbeta_ds[0]+factor2*dP_ds[0]+factor3*da_ds[0];
						dF_dsigma[1] = factor1*dbeta_ds[1]+factor2*dP_ds[1]+factor3*da_ds[1];
						dF_dsigma[2] = factor1*dbeta_ds[2]+factor2*dP_ds[2]+factor3*da_ds[2];
						dF_dsigma[3] = factor1*dbeta_ds[3]+factor2*dP_ds[3]+factor3*da_ds[3];
						dF_dsigma[4] = factor1*dbeta_ds[4]+factor2*dP_ds[4]+factor3*da_ds[4];
						dF_dsigma[5] = factor1*dbeta_ds[5]+factor2*dP_ds[5]+factor3*da_ds[5];

						if (rd2F_d2Sigma!=0)
						{
							for (int j=0; j<6;j++)
							{
								for (int i=j; i<6;i++)
								{
									double tmp = (factor5*dbeta_ds[j]+factor4*dP_ds[j]) * dbeta_ds[i] + factor1*dbeta_ds2[6*j+i]
												  +factor4*dbeta_ds[j]*dP_ds[i] + factor2*dP_ds2[6*j+i];
									d2F_d2Sigma[6*j+i] = tmp;
									if (i!=j)
									{
										d2F_d2Sigma[6*i+j] = tmp;
									}
								}
							}
						}
					}
					return F;
				}
			}// D<=-tolerance_D
			else
			{
				// two identical solutions with principal[1] and principal[2], since only one is positive, the solution is principal[0]
				F = principal[0]-rFct;
				if (rdF_dSigma!=0)
				{
					double df_dq,
					df_da,
					df_dP,
					df_dq2,
					df_dP2,
					df_dqdP,


					P2 = P*P;
					P3 = P2*P;

					df_dq  = -2./(9.*P2);
					df_dP  = -4./3.;
					df_da  = -1./3.;

					df_dq2  = 16./(243*P2*P3);
					df_dqdP = 20./(81*P3);;
					df_dP2  = -20./(27.*P);

					dF_dsigma[0] = df_dq*dq_ds[0]+df_dP*dP_ds[0]+df_da*da_ds[0];
					dF_dsigma[1] = df_dq*dq_ds[1]+df_dP*dP_ds[1]+df_da*da_ds[1];
					dF_dsigma[2] = df_dq*dq_ds[2]+df_dP*dP_ds[2]+df_da*da_ds[2];
					dF_dsigma[3] = df_dq*dq_ds[3]+df_dP*dP_ds[3]+df_da*da_ds[3];
					dF_dsigma[4] = df_dq*dq_ds[4]+df_dP*dP_ds[4]+df_da*da_ds[4];
					dF_dsigma[5] = df_dq*dq_ds[5]+df_dP*dP_ds[5]+df_da*da_ds[5];

					if (rd2F_d2Sigma!=0)
					{
						for (int j=0; j<6;j++)
						{
							for (int i=j; i<6;i++)
							{
								double tmp = (df_dq2*dq_ds[j]+df_dqdP*dP_ds[j])*dq_ds[i]+df_dq*dq_ds2[6*j+i]
													 +(df_dqdP*dq_ds[j]+df_dP2*dP_ds[j])*dP_ds[i]+df_dP*dP_ds2[6*j+i];
								d2F_d2Sigma[6*j+i] = tmp;
								if (i!=j)
								{
									d2F_d2Sigma[6*i+j] = tmp;
								}
							}
						}
					}
				}
				return F;
			}//else D<=-tolerance_D
		}// num_pos==1
		else
		{
			assert(num_pos==2);
			if (D<=-tolerance_D)
			{
				double dprinc_dsigma[3][6],d2princ_d2sigma[3][36];
				F = 0;
				if (principal[0]>0)
				{
					F += principal[0]*principal[0];
					if (rdF_dSigma!=0)
					{
						// calculate derivatives
						sin_beta = sin(beta);
						factor4 = 2.*sin_beta;
						factor1 = P*factor4;
						factor2 = -2.*cos_beta;
						factor3 = -1./3.;
						factor5 = -P*factor2;
						dprinc_dsigma[0][0] = factor1*dbeta_ds[0]+factor2*dP_ds[0]+factor3*da_ds[0];
						dprinc_dsigma[0][1] = factor1*dbeta_ds[1]+factor2*dP_ds[1]+factor3*da_ds[1];
						dprinc_dsigma[0][2] = factor1*dbeta_ds[2]+factor2*dP_ds[2]+factor3*da_ds[2];
						dprinc_dsigma[0][3] = factor1*dbeta_ds[3]+factor2*dP_ds[3]+factor3*da_ds[3];
						dprinc_dsigma[0][4] = factor1*dbeta_ds[4]+factor2*dP_ds[4]+factor3*da_ds[4];
						dprinc_dsigma[0][5] = factor1*dbeta_ds[5]+factor2*dP_ds[5]+factor3*da_ds[5];

						if (rd2F_d2Sigma!=0)
						{
							for (int j=0; j<6;j++)
							{
								for (int i=j; i<6;i++)
								{
									double tmp = (factor5*dbeta_ds[j]+factor4*dP_ds[j]) * dbeta_ds[i] + factor1*dbeta_ds2[6*j+i]
																+factor4*dbeta_ds[j]*dP_ds[i] + factor2*dP_ds2[6*j+i];
									d2princ_d2sigma[0][6*j+i] = tmp;
									if (i!=j)
									{
										d2princ_d2sigma[0][6*i+j] = tmp;
									}
								}
							}
						}
					}
				}
				if (principal[1]>0)
				{
					F += principal[1]*principal[1];
					if (rdF_dSigma!=0)
					{
						sin_beta_add = sin(beta+M_PI/3.);
						factor4 = -2.*sin_beta_add;
						factor1 = P*factor4;
						factor2 = 2.*cos_beta_add;
						factor3 = -1./3.;
						factor5 = -P*factor2;
						dprinc_dsigma[1][0] = factor1*dbeta_ds[0]+factor2*dP_ds[0]+factor3*da_ds[0];
						dprinc_dsigma[1][1] = factor1*dbeta_ds[1]+factor2*dP_ds[1]+factor3*da_ds[1];
						dprinc_dsigma[1][2] = factor1*dbeta_ds[2]+factor2*dP_ds[2]+factor3*da_ds[2];
						dprinc_dsigma[1][3] = factor1*dbeta_ds[3]+factor2*dP_ds[3]+factor3*da_ds[3];
						dprinc_dsigma[1][4] = factor1*dbeta_ds[4]+factor2*dP_ds[4]+factor3*da_ds[4];
						dprinc_dsigma[1][5] = factor1*dbeta_ds[5]+factor2*dP_ds[5]+factor3*da_ds[5];

						if (rd2F_d2Sigma!=0)
						{
							for (int j=0; j<6;j++)
							{
								for (int i=j; i<6;i++)
								{
									double tmp = (factor5*dbeta_ds[j]+factor4*dP_ds[j]) * dbeta_ds[i] + factor1*dbeta_ds2[6*j+i]
																+factor4*dbeta_ds[j]*dP_ds[i] + factor2*dP_ds2[6*j+i];
									d2princ_d2sigma[1][6*j+i] = tmp;
									if (i!=j)
									{
										d2princ_d2sigma[1][6*i+j] = tmp;
									}
								}
							}
						}
					}
				}
				if (principal[2]>0)
				{
					F += principal[2]*principal[2];
					if (rdF_dSigma!=0)
					{
						sin_beta_sub = sin(beta-M_PI/3.);
						factor4 = -2.*sin_beta_sub;
						factor1 = P*factor4;
						factor2 = 2.*cos_beta_sub;
						factor3 = -1./3.;
						factor5 = -P*factor2;
						dprinc_dsigma[2][0] = factor1*dbeta_ds[0]+factor2*dP_ds[0]+factor3*da_ds[0];
						dprinc_dsigma[2][1] = factor1*dbeta_ds[1]+factor2*dP_ds[1]+factor3*da_ds[1];
						dprinc_dsigma[2][2] = factor1*dbeta_ds[2]+factor2*dP_ds[2]+factor3*da_ds[2];
						dprinc_dsigma[2][3] = factor1*dbeta_ds[3]+factor2*dP_ds[3]+factor3*da_ds[3];
						dprinc_dsigma[2][4] = factor1*dbeta_ds[4]+factor2*dP_ds[4]+factor3*da_ds[4];
						dprinc_dsigma[2][5] = factor1*dbeta_ds[5]+factor2*dP_ds[5]+factor3*da_ds[5];

						if (rd2F_d2Sigma!=0)
						{
							for (int j=0; j<6;j++)
							{
								for (int i=j; i<6;i++)
								{
									double tmp = (factor5*dbeta_ds[j]+factor4*dP_ds[j]) * dbeta_ds[i] + factor1*dbeta_ds2[6*j+i]
																+factor4*dbeta_ds[j]*dP_ds[i] + factor2*dP_ds2[6*j+i];
									d2princ_d2sigma[2][6*j+i] = tmp;
									if (i!=j)
									{
										d2princ_d2sigma[2][6*i+j] = tmp;
									}
								}
							}
						}
					}
				}
				int s1,s2;
				double f;
				// calculate yield function
				f = sqrt(F);
				F = f-rFct;
				if (rdF_dSigma!=0)
				{
					if (principal[0]>0)
					{
						s1 = 0;
						if (principal[0]>0)
						{
							s2 = 1;
						}
						else
						{
							s2 = 2;
						}
					}
					else
					{
						s1 = 1;
						s2 = 2;
					}

					// calculate gradient and Hessian
					factor1 = 1./f;
					dF_dsigma[0] = factor1*(dprinc_dsigma[s1][0]*principal[s1] + dprinc_dsigma[s2][0]*principal[s2]);
					dF_dsigma[1] = factor1*(dprinc_dsigma[s1][1]*principal[s1] + dprinc_dsigma[s2][1]*principal[s2]);
					dF_dsigma[2] = factor1*(dprinc_dsigma[s1][2]*principal[s1] + dprinc_dsigma[s2][2]*principal[s2]);
					dF_dsigma[3] = factor1*(dprinc_dsigma[s1][3]*principal[s1] + dprinc_dsigma[s2][3]*principal[s2]);
					dF_dsigma[4] = factor1*(dprinc_dsigma[s1][4]*principal[s1] + dprinc_dsigma[s2][4]*principal[s2]);
					dF_dsigma[5] = factor1*(dprinc_dsigma[s1][5]*principal[s1] + dprinc_dsigma[s2][5]*principal[s2]);

					if (rd2F_d2Sigma!=0)
					{
						factor2 = factor1*factor1*factor1;
						factor3 = (f*f);

						factor1 = principal[s1]*principal[s1];
						factor4 = principal[s2]*principal[s2];
						factor5 = principal[s1]*principal[s2];
						for (int j=0; j<6;j++)
						{
							for (int i=j; i<6;i++)
							{
								double tmp = factor2*(
														 d2princ_d2sigma[s1][6*j+i]*factor3*principal[s1]
														 +d2princ_d2sigma[s2][6*j+i]*factor3*principal[s2]
														 +dprinc_dsigma[s1][i]*(factor4*dprinc_dsigma[s1][j]-factor5*dprinc_dsigma[s2][j])
														 +dprinc_dsigma[s2][i]*(factor1*dprinc_dsigma[s2][j]-factor5*dprinc_dsigma[s1][j])
													 );
								d2F_d2Sigma[6*j+i] = tmp;
								if (i!=j)
								{
									d2F_d2Sigma[6*i+j] = tmp;
								}
							}
						}
					}
				}
				return F;

			}// D<=-tolerance_D
			else
			{
				// two identical solutions with principal[1] and principal[2]
				//F = sqrt(2.)*(principal[1])-rFct; this is not 100% correct, since they are only equal up to a certain accuracy
				F = sqrt(principal[1]*principal[1]+principal[2]*principal[2])-rFct;

				if (rdF_dSigma!=0)
				{
					double df_dq,
					df_da,
					df_dP,
					df_dq2,
					df_dP2,
					df_dqdP,
					df_dqda,
					df_dPda,
					mul_3P_sub_a,
					sign_3P_sub_a,
					mul_3P_sub_a_pow_2,
					mul_3P_sub_a_pow_3;

					mul_3P_sub_a = 3.*P-a;
					mul_3P_sub_a_pow_2 = mul_3P_sub_a*mul_3P_sub_a;
					mul_3P_sub_a_pow_3 = mul_3P_sub_a_pow_2*mul_3P_sub_a;
					sign_3P_sub_a = mul_3P_sub_a<0 ? -1 : 1;

					P2 = P*P;
					P3 = P2*P;
					P6 = P3*P3;

					df_dq  = -sqrt(2)*(6.*P+a)*sign_3P_sub_a/(9.*P2*mul_3P_sub_a);
					df_dP  = sqrt(2)*(15.*P-2.*a)*sign_3P_sub_a/(3.*mul_3P_sub_a);
					df_da  = -sqrt(2)/(3.)*sign_3P_sub_a;

					df_dq2  = -sqrt(2)*(54.*P3+216.*P2*a+27.*a2*P-8.*a3)*sign_3P_sub_a/(243.*P2*P3*mul_3P_sub_a_pow_3);
					df_dqdP = sqrt(2)*(1026.*P3-27.*P2*a-54.*a2*P+10.*a3)*sign_3P_sub_a/(81.*P3*mul_3P_sub_a_pow_3);
					df_dqda = -sqrt(2)*sign_3P_sub_a/(mul_3P_sub_a_pow_2*P);
					df_dP2  = -sqrt(2)*(1026.*P3+216.*P2*a-135.*a2*P+10.*a3)*sign_3P_sub_a/(27.*P*mul_3P_sub_a_pow_3);
					df_dPda = sqrt(2)*(3.*P)*sign_3P_sub_a/(mul_3P_sub_a_pow_2);

					dF_dsigma[0] = df_dq*dq_ds[0]+df_dP*dP_ds[0]+df_da*da_ds[0];
					dF_dsigma[1] = df_dq*dq_ds[1]+df_dP*dP_ds[1]+df_da*da_ds[1];
					dF_dsigma[2] = df_dq*dq_ds[2]+df_dP*dP_ds[2]+df_da*da_ds[2];
					dF_dsigma[3] = df_dq*dq_ds[3]+df_dP*dP_ds[3]+df_da*da_ds[3];
					dF_dsigma[4] = df_dq*dq_ds[4]+df_dP*dP_ds[4]+df_da*da_ds[4];
					dF_dsigma[5] = df_dq*dq_ds[5]+df_dP*dP_ds[5]+df_da*da_ds[5];

					if (rd2F_d2Sigma!=0)
					{
						for (int j=0; j<6;j++)
						{
							for (int i=j; i<6;i++)
							{
								double tmp = (df_dq2*dq_ds[j]+df_dqdP*dP_ds[j]+df_dqda*da_ds[j])*dq_ds[i]+df_dq*dq_ds2[6*j+i]
													 +(df_dqdP*dq_ds[j]+df_dP2*dP_ds[j]+df_dPda*da_ds[j])*dP_ds[i]+df_dP*dP_ds2[6*j+i]
													 +(df_dqda*dq_ds[j]+df_dPda*dP_ds[j])*da_ds[i];
								d2F_d2Sigma[6*j+i] = tmp;
								if (i!=j)
								{
									d2F_d2Sigma[6*i+j] = tmp;

								}
							}
						}
					}
				}
				return F;
			}
		}// num_pos==2
		break; //num_pos==1,2
	case 3:
	{
		double f = sqrt(principal[0]*principal[0]+principal[1]*principal[1]+principal[2]*principal[2]);
		F = f-rFct;
		if (rdF_dSigma!=0)
		{
			// a
			da_ds[0] = -1.;
			da_ds[1] = -1.;
			da_ds[2] = -1.;
			da_ds[3] =  0.;
			da_ds[4] =  0.;
			da_ds[5] =  0.;

			// b
			db_ds[0] = rStress(1)+rStress(2);
			db_ds[1] = rStress(0)+rStress(2);
			db_ds[2] = rStress(0)+rStress(1);
			db_ds[3] = -2.*rStress(3);
			db_ds[4] = -2.*rStress(4);
			db_ds[5] = -2.*rStress(5);

			db_ds2[1]  = 1.;
			db_ds2[2]  = 1.;
			db_ds2[8]  = 1.;
			db_ds2[21]  = -2.;
			db_ds2[28] = -2.;
			db_ds2[35] = -2.;

			double df_da,df_db,df_da2,df_db2,df_dadb;

			factor1 = 1./f;
			factor3 = factor1/(f*f);

			df_da = a*factor1;
			df_db = -factor1;

			df_da2 = -2.*b*factor3;
			df_dadb = a*factor3;
			df_db2 = -factor3;

			dF_dsigma[0] = df_da*da_ds[0]+df_db*db_ds[0];
			dF_dsigma[1] = df_da*da_ds[1]+df_db*db_ds[1];
			dF_dsigma[2] = df_da*da_ds[2]+df_db*db_ds[2];
			dF_dsigma[3] = df_da*da_ds[3]+df_db*db_ds[3];
			dF_dsigma[4] = df_da*da_ds[4]+df_db*db_ds[4];
			dF_dsigma[5] = df_da*da_ds[5]+df_db*db_ds[5];

			if (rd2F_d2Sigma!=0)
			{
				for (int j=0; j<6;j++)
				{
					for (int i=j; i<6;i++)
					{
						double tmp = (df_da2*da_ds[j]+df_dadb*db_ds[j])*da_ds[i]
											 +(df_dadb*da_ds[j]+df_db2*db_ds[j])*db_ds[i]+df_db*db_ds2[6*j+i];
						d2F_d2Sigma[6*j+i] = tmp;
						if (i!=j)
						{
							d2F_d2Sigma[6*i+j] = tmp;
						}
					}
				}
			}
		}
		return F;
	}
	break;
	default:
		throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine3DRounded] programming error.");
	}
	throw MechanicsException("[NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceRankine3DRounded] programming error (end of routine).");
}

//! @brief calculates the Drucker Prager yield surface and the derivatives with respect to the stress
//! @param rStress current stress
//! @param rBETA parameter of the Drucker Prager yield surface
//! @param rHP parameter of the Drucker Prager yield surface
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rErrorDerivatives true, if derivative can't be calculated (on the hydrostatic axis)
//! @return yield function
double NuTo::GradientDamagePlasticityEngineeringStress::YieldSurfaceDruckerPrager3D(Eigen::Matrix<double,6,1>& rStress, double rBeta, double rHP,
		Eigen::Matrix<double,6,1>* rdF_dSigma,Eigen::Matrix<double,6,6>* rd2F_d2Sigma, bool &rErrorDerivatives
        )const
{
    double stress3_pow2 = rStress(3)*rStress(3);
    double stress4_pow2 = rStress(4)*rStress(4);
    double stress5_pow2 = rStress(5)*rStress(5);

    // *******************************************************************
    // 	first invariante							 					 *
    // *******************************************************************
    double invariante_1 = rStress(0)+rStress(1)+rStress(2);


    // *******************************************************************
    // 	second invariante							 					 *
    // *******************************************************************
    double invariante_2 = ((rStress(0)-rStress(1))*(rStress(0)-rStress(1))+
                    (rStress(1)-rStress(2))*(rStress(1)-rStress(2))+
                    (rStress(0)-rStress(2))*(rStress(0)-rStress(2)))/6.+
                   rStress(3)*rStress(3)+
                   rStress(4)*rStress(4)+
                   rStress(5)*rStress(5);

    // *******************************************************************
    // 	F_DEV:	term of stress-deviator with/without kinematic hard.     *
    // *******************************************************************
    double F_DEV = sqrt( invariante_2 );


    // *******************************************************************
    // 	yield value     *
    // *******************************************************************
    double F_DP  = invariante_1 * rBeta/3. + F_DEV - rHP ;

    if (rdF_dSigma!=0)
    {
        if (invariante_2>=1e-10)
        {
            if (rdF_dSigma!=0)
            {
            	double *dF_dsigma =  (*rdF_dSigma).data();

				// gradient
				double factor = 1./(F_DEV*6.);
				dF_dsigma[0] = factor * (2.*rStress(0)-rStress(1)-rStress(2))+rBeta/3.;
				dF_dsigma[1] = factor * (2.*rStress(1)-rStress(0)-rStress(2))+rBeta/3.;
				dF_dsigma[2] = factor * (2.*rStress(2)-rStress(0)-rStress(1))+rBeta/3.;
				dF_dsigma[3] = factor * (6.*rStress(3)); /* vector notation from second order tensor */
				dF_dsigma[4] = factor * (6.*rStress(4)); /* vector notation from second order tensor */
				dF_dsigma[5] = factor * (6.*rStress(5)); /* vector notation from second order tensor */
            }

			if (rd2F_d2Sigma!=0)
			{
				double *d2F_d2sigma =  (*rd2F_d2Sigma).data();
				//hessian
				double factor = 1./(invariante_2*F_DEV*12.);
				double tmp_scalar = rStress(1)-rStress(2);
				d2F_d2sigma[0] = factor * (tmp_scalar*tmp_scalar+4.*(stress3_pow2+stress4_pow2+stress5_pow2));
				d2F_d2sigma[1] = factor * (rStress(2)*(rStress(0)+rStress(1)-rStress(2))-rStress(0)*rStress(1)-2.*(stress3_pow2+stress4_pow2+stress5_pow2));
				d2F_d2sigma[2] = factor * (rStress(1)*(rStress(0)+rStress(2)-rStress(1))-rStress(0)*rStress(2)-2.*(stress3_pow2+stress4_pow2+stress5_pow2));
				tmp_scalar = -2.*(2*rStress(0)-rStress(1)-rStress(2));
				d2F_d2sigma[3] = factor * rStress(3)*tmp_scalar;
				d2F_d2sigma[4] = factor * rStress(4)*tmp_scalar;
				d2F_d2sigma[5] = factor * rStress(5)*tmp_scalar;

				d2F_d2sigma[6] = d2F_d2sigma[1];
				tmp_scalar = rStress(0)-rStress(2);
				d2F_d2sigma[7] = factor * (tmp_scalar*tmp_scalar+4.*(stress3_pow2+stress4_pow2+stress5_pow2));
				d2F_d2sigma[8] = factor * (rStress(0)*(rStress(1)+rStress(2)-rStress(0))-rStress(1)*rStress(2)-2.*(stress3_pow2+stress4_pow2+stress5_pow2));
				tmp_scalar = 2.*(rStress(0)+rStress(2)-2*rStress(1));
				d2F_d2sigma[9] = factor * rStress(3)*tmp_scalar;
				d2F_d2sigma[10] = factor * rStress(4)*tmp_scalar;
				d2F_d2sigma[11] = factor * rStress(5)*tmp_scalar;

				d2F_d2sigma[12] = d2F_d2sigma[2];
				d2F_d2sigma[13] = d2F_d2sigma[8];
				tmp_scalar = rStress(0)-rStress(1);
				d2F_d2sigma[14] = factor * (tmp_scalar*tmp_scalar+4.*(stress3_pow2+stress4_pow2+stress5_pow2));
				tmp_scalar = 2.*(rStress(0)+rStress(1)-2*rStress(2));
				d2F_d2sigma[15] = factor * rStress(3)*tmp_scalar;
				d2F_d2sigma[16] = factor * rStress(4)*tmp_scalar;
				d2F_d2sigma[17] = factor * rStress(5)*tmp_scalar;

				tmp_scalar = 4.*(rStress(0)*rStress(0) + rStress(1)*rStress(1) + rStress(2)*rStress(2)
								 - rStress(0)*rStress(1) - rStress(0)*rStress(2) - rStress(1)*rStress(2));
				d2F_d2sigma[18] = d2F_d2sigma[3];
				d2F_d2sigma[19] = d2F_d2sigma[9];
				d2F_d2sigma[20] = d2F_d2sigma[15];
				d2F_d2sigma[21] = factor*(tmp_scalar+12.*(stress4_pow2+stress5_pow2));
				d2F_d2sigma[22] = factor*(-12.*rStress(3)*rStress(4));
				d2F_d2sigma[23] = factor*(-12.*rStress(3)*rStress(5));

				d2F_d2sigma[24] = d2F_d2sigma[4];
				d2F_d2sigma[25] = d2F_d2sigma[10];
				d2F_d2sigma[26] = d2F_d2sigma[16];
				d2F_d2sigma[27] = d2F_d2sigma[22];
				d2F_d2sigma[28] = factor*(tmp_scalar+12.*(stress3_pow2+stress5_pow2));
				d2F_d2sigma[29] = factor*(-12.*rStress(4)*rStress(5));

				d2F_d2sigma[30] = d2F_d2sigma[5];
				d2F_d2sigma[31] = d2F_d2sigma[11];
				d2F_d2sigma[32] = d2F_d2sigma[17];
				d2F_d2sigma[33] = d2F_d2sigma[23];
				d2F_d2sigma[34] = d2F_d2sigma[29];
				d2F_d2sigma[35] = factor*(tmp_scalar+12.*(stress3_pow2+stress4_pow2));
            }
        }
        else
        {
            rErrorDerivatives = false;
        }
    }
    return F_DP;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
//! @param rLogger stream for the output
void NuTo::GradientDamagePlasticityEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus      : " << this->mE << "\n";
    rLogger << "    Poisson's ratio      : " << this->mNu << "\n";
    rLogger << "    nonlocal radius      : " << this->mNonlocalRadius << "\n";
    rLogger << "    tensile strength     : " << this->mTensileStrength << "\n";
    rLogger << "    compressive strength : " << this->mCompressiveStrength << "\n";
    rLogger << "    biaxial compressive strength : " << this->mBiaxialCompressiveStrength << "\n";
    rLogger << "    fracture energy      : " << this->mFractureEnergy << "\n";
    rLogger << "    thermal expansion coeff : " << this->mThermalExpansionCoefficient << "\n";
}

// check parameters
void NuTo::GradientDamagePlasticityEngineeringStress::CheckParameters()const
{
    this->CheckBiaxialCompressiveStrength(this->mBiaxialCompressiveStrength);
    this->CheckDensity(this->mRho);
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckNonlocalRadius(this->mNonlocalRadius);
    this->CheckTensileStrength(this->mTensileStrength);
    this->CheckCompressiveStrength(this->mCompressiveStrength);
    this->CheckFractureEnergy(this->mFractureEnergy);
    this->CheckThermalExpansionCoefficient(this->mThermalExpansionCoefficient);
}

double NuTo::GradientDamagePlasticityEngineeringStress::CalculateDamage(double rkappa, double rKappaD)const
{
    double omega = 1.-exp(-rkappa/rKappaD);
    if (omega>MAX_OMEGA)
    {
        omega = MAX_OMEGA;
    }
    return omega;
}

double NuTo::GradientDamagePlasticityEngineeringStress::CalculateDerivativeDamage(double rkappa, double rKappaD, double& rdOmegadKappa)const
{
    double invKappaD(1./rKappaD);
    double tmp = exp(-rkappa*invKappaD);
    if (1.-tmp>MAX_OMEGA)
    {
        rdOmegadKappa = 0.;
        return MAX_OMEGA;
    }
    else
    {
        rdOmegadKappa = invKappaD * tmp;
        return 1.-tmp;
    }
}


double NuTo::GradientDamagePlasticityEngineeringStress::CalculateKappaD()const
{
    return 15.*mFractureEnergy/(16.*mTensileStrength);
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::GradientDamagePlasticityEngineeringStress::HaveTmpStaticData() const
{
    return false;
}

//! @brief ... returns true, if a material model has is nonlocal (stiffness is of dynamic size, nonlocal averaging)
//! @return ... see brief explanation
bool NuTo::GradientDamagePlasticityEngineeringStress::IsNonlocalModel()const
{
    if (mDamage)
        return true;
    else
        return false;
}
