// $Id$
// NonlocalDamagePlasticityEngineeringStress.cpp
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

#include "base/Logger.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.h"
#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity3D.h"
#include "mechanics/constitutive/mechanics/Damage.h"
#include "mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "mechanics/constitutive/mechanics/NonlocalDamagePlasticityEngineeringStress.h"
#include "mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress1D.h"
#include "mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress2D.h"
#include "mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress3D.h"
#include "mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/sections/SectionBase.h"
#include "mechanics/sections/SectionEnum.h"

#define sqrt3 1.732050808
#define MAX_OMEGA 0.999
//#define ENABLE_DEBUG

NuTo::NonlocalDamagePlasticityEngineeringStress::NonlocalDamagePlasticityEngineeringStress() : ConstitutiveBase()
{
    mRho = 0.;
    mE = 0.;
    mNu = 0.;
    mNonlocalRadius = 1.;
    mTensileStrength = 0.;
    mCompressiveStrength = 0.;
    mBiaxialCompressiveStrength = 0.;
    mFractureEnergy = 0.;
    mYieldSurface = Constitutive::COMBINED_ROUNDED;
    mDamage = true;
    mThermalExpansionCoefficient = 0.;
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NonlocalDamagePlasticityEngineeringStress::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
       std::cout << "start serialize NonlocalDamagePlasticityEngineeringStress" << "\n";
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
       std::cout << "finish serialize NonlocalDamagePlasticityEngineeringStress" << "\n";
#endif
    }
    BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NonlocalDamagePlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate1D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate1D] not implemented for 1D.");
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate2D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    // get interpolation type information determining which input on the constitutive level should be used
    const InterpolationType* interpolationType = rElement->GetInterpolationType();

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }

    if (rElement->GetSection()->GetType()!=Section::PLANE_STRAIN)
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate] only plane strain is implemented.");

    EngineeringStrain2D engineeringStrain;
    // calculate engineering strain
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D)==rConstitutiveInput.end())
        throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate] deformation gradient 2d needed to evaluate engineering strain2d.");
    const DeformationGradient2D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D)->second->GetDeformationGradient2D());
    deformationGradient.GetEngineeringStrain(engineeringStrain);

    //a recalculation of the plastic strain is not necessary, since this has been performed at the previous iteration with the update of the nonlocaltmpstatic data
    //Get previous ip_data
    ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain *oldStaticData = (rElement->GetStaticData(rIp))->AsNonlocalDamagePlasticity2DPlaneStrain();

    // subtract local plastic strain that has been calculated within the updatetmpstaticdata routine
    double elastStrain[4];
    elastStrain[0] = engineeringStrain[0] - oldStaticData->mTmpEpsilonP[0];
    elastStrain[1] = engineeringStrain[1] - oldStaticData->mTmpEpsilonP[1];
    elastStrain[2] = engineeringStrain[2] - oldStaticData->mTmpEpsilonP[2];
    elastStrain[3] =  - oldStaticData->mTmpEpsilonP[3];

    // subtract thermal strain
    if (interpolationType->IsConstitutiveInput(Node::TEMPERATURE))
    {
        std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>::const_iterator itInput(rConstitutiveInput.find(NuTo::Constitutive::Input::TEMPERATURE));
        if (itInput==rConstitutiveInput.end())
            throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate2D] temperature needed to evaluate thermal engineering strain2d.");
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
        engineeringStress[0] = oneMinusOmega*(C11 * elastStrain[0] + C12*(elastStrain[1]+elastStrain[3]));
        engineeringStress[1] = oneMinusOmega*(C11 * elastStrain[1] + C12*(elastStrain[0]+elastStrain[3]));
        engineeringStress[2] = oneMinusOmega*(C33*elastStrain[2]);
    }
    else
    {
    	engineeringStress[0] = (C11 * elastStrain[0] + C12*(elastStrain[1]+elastStrain[3]));
    	engineeringStress[1] = (C11 * elastStrain[1] + C12*(elastStrain[0]+elastStrain[3]));
    	engineeringStress[2] = (C33*elastStrain[2]);
    }

    //set this to true, if update is in the map, perform the update after all other outputs have been calculated
    bool performUpdateAtEnd(false);

    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
            itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch(itOutput->first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_2D:
        {
            // copy Engineering stress
            itOutput->second->GetEngineeringStress2D() = engineeringStress;

        }
        break;
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
        {
            //this is for the visualize routines
            EngineeringStress3D& engineeringStress3D(itOutput->second->GetEngineeringStress3D());

            engineeringStress3D[0] = engineeringStress[0];
            engineeringStress3D[1] = engineeringStress[1];
            engineeringStress3D[2] = oneMinusOmega*(C11 * elastStrain[3] + C12*(elastStrain[0]+elastStrain[1]));
            engineeringStress3D[3] = engineeringStress[2];
            engineeringStress3D[4] = 0.;
            engineeringStress3D[5] = 0.;
        }
        break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D:
        {
            if (mDamage)
            {
                //calculate damage parameter from the equivalente plastic strain
                //it is scaled related to the length (based on the direction of the first principal plastic strain)
                double dOmegadKappa;
                double oneMinusOmega = 1.-CalculateDerivativeDamage(kappa, kappaD, dOmegadKappa);
                assert(std::abs(oneMinusOmega+omega-1.)<1e-10);
                //rLogger << "omega " << 1.-oneMinusOmega << "\n";
                if (unloading)
                {
                	ConstitutiveOutputBase* outputBase(itOutput->second);
                	outputBase->SetLocalSolution(true);
                	ConstitutiveTangentLocal<3,3>& tangent(outputBase->GetSubMatrix_3x3(0).AsConstitutiveTangentLocal_3x3());
                    //rLogger << "unloading local stiffness " << C11 << " "<< C12 << " " << C33 << " omega " << 1.-oneMinusOmega << "\n";

                    tangent.SetSymmetry(true);
                    tangent(0,0) = oneMinusOmega * C11;
                    tangent(0,1) = oneMinusOmega * C12;
                    tangent(0,2) = 0.;

                    tangent(1,0) = tangent(0,1);
                    tangent(1,1) = tangent(0,0);
                    tangent(1,2) = 0.;

                    tangent(2,0) = 0.;
                    tangent(2,1) = 0.;
                    tangent(2,2) = oneMinusOmega * C33;
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

                            tangent = sigmaElast * minusdOmegadEpsilon.transpose() ;

                            if (nonlocalElements[countNonlocalElement]==rElement && rIp == theNonlocalIP)
                            {

                                const double *tmpdEP(oldStaticData->mTmpdEpsilonPdEpsilon);
                                tangent(0,0) += oneMinusOmega*(C11*(1.-tmpdEP[0]) - C12*(tmpdEP[1] + tmpdEP[3]));
                                tangent(1,0) += oneMinusOmega*(C12*(1.-tmpdEP[0] - tmpdEP[3]) - C11*tmpdEP[1]);
                                tangent(2,0) += oneMinusOmega*(-C33*tmpdEP[2]);
                                tangent(0,1) += oneMinusOmega*(-C11*tmpdEP[4] + C12*(1.-tmpdEP[5]- tmpdEP[7]));
                                tangent(1,1) += oneMinusOmega*(-C12*(tmpdEP[4] + tmpdEP[7]) + C11*(1.-tmpdEP[5]));
                                tangent(2,1) += oneMinusOmega*(-C33*tmpdEP[6]);
                                tangent(0,2) += oneMinusOmega*(-C11*tmpdEP[8] - C12*(tmpdEP[9]+tmpdEP[11]));
                                tangent(1,2) += oneMinusOmega*(-C12*(tmpdEP[8]+ tmpdEP[11]) - C11*tmpdEP[9]);
                                tangent(2,2) += oneMinusOmega*(C33*(1.-tmpdEP[10]));
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

                const double *tmpdEP(oldStaticData->mTmpdEpsilonPdEpsilon);
                tangent(0,0) = C11*(1.-tmpdEP[0]) - C12*(tmpdEP[1] + tmpdEP[3]);
                tangent(1,0) = C12*(1.-tmpdEP[0] - tmpdEP[3]) - C11*tmpdEP[1];
                tangent(2,0) = -C33*tmpdEP[2];
                tangent(0,1) = -C11*tmpdEP[4] + C12*(1.-tmpdEP[5]- tmpdEP[7]);
                tangent(1,1) = -C12*(tmpdEP[4] + tmpdEP[7]) + C11*(1.-tmpdEP[5]);
                tangent(2,1) = -C33*tmpdEP[6];
                tangent(0,2) = -C11*tmpdEP[8] - C12*(tmpdEP[9]+tmpdEP[11]);
                tangent(1,2) = -C12*(tmpdEP[8]+ tmpdEP[11]) - C11*tmpdEP[9];
                tangent(2,2) = C33*(1.-tmpdEP[10]);
            }
        }
        break;
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
        {
            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
			engineeringStrain3D[0] = engineeringStrain[0];
			engineeringStrain3D[1] = engineeringStrain[1];
			engineeringStrain3D[2] = 0;
			engineeringStrain3D[3] = engineeringStrain[2];
			engineeringStrain3D[4] = 0.;
			engineeringStrain3D[5] = 0.;
        }
        break;
        case NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D:
        {
        	EngineeringStrain3D& engineeringPlasticStrain(itOutput->second->GetEngineeringStrain3D());
        	engineeringPlasticStrain[0] = oldStaticData->mTmpEpsilonP[0];
        	engineeringPlasticStrain[1] = oldStaticData->mTmpEpsilonP[1];
        	engineeringPlasticStrain[2] = oldStaticData->mTmpEpsilonP[3];
        	engineeringPlasticStrain[3] = oldStaticData->mTmpEpsilonP[2];
        	engineeringPlasticStrain[4] = 0.;
        	engineeringPlasticStrain[5] = 0.;
        }
        break;
        case NuTo::Constitutive::Output::DAMAGE:
        {
            itOutput->second->GetDamage().SetDamage(omega);        }
        break;
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
        	if (rConstitutiveOutput.size()!=1)
        		throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate3D] tmp_static_data has to be updated without any other outputs, call it separately.");

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

			//std::cout << "new plastic strain " << rNewEpsilonP << std::endl;

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
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
        	performUpdateAtEnd = true;
        }
        break;
        default:
            throw MechanicsException(std::string("[NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate3D] output object)") +
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
            (engineeringStress[0]+prevStress[0])*(engineeringStrain[0]-prevStrain[0])+
            (engineeringStress[1]+prevStress[1])*(engineeringStrain[1]-prevStrain[1])+
            (engineeringStress[2]+prevStress[2])*(engineeringStrain[2]-prevStrain[2])
            );
        oldStaticData->mPrevStrain = engineeringStrain;
        oldStaticData->mPrevSigma = engineeringStress;
        oldStaticData->SetPrevTotalEnergy(energy);

        // update the parts of the static data that are not related to the temporary updates (from the nonlocal calculation)
        oldStaticData->mKappa = oldStaticData->mTmpKappa;

//        Eigen::Matrix<double,4,1>::Map(oldStaticData->mEpsilonP,4,1) = Eigen::Matrix<double,4,1>::Map(oldStaticData->mTmpEpsilonP,4,1);;
    }

    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate3D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::Evaluate3D] not implemented for 3D.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
        const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for 1D.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
        const ElementBase* rElement) const
{
    if (rElement->GetSection()==0)
        throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Section required to distinguish between plane stress and plane strain and thickness information.");
    if (rElement->GetSection()->GetType()==NuTo::Section::PLANE_STRESS)
        throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for plane stress.");
    else
        return new ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain();
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticityEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
        const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataNonlocalDamagePlasticity3D();
}


// calculate coefficients of the material matrix
void NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE/((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE/(2*(1.0 + this->mNu));
}

// parameters /////////////////////////////////////////////////////////////

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::NonlocalDamagePlasticityEngineeringStress::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        return mCompressiveStrength;
    case Constitutive::eConstitutiveParameter::BIAXIAL_COMPRESSIVE_STRENGTH:
        return mBiaxialCompressiveStrength;
    case Constitutive::eConstitutiveParameter::DENSITY:
        return this->mRho;
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        return mFractureEnergy;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        return this->mNonlocalRadius;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        return this->mNu;
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
        return mTensileStrength;
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        return this->mThermalExpansionCoefficient;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        return this->mE;
    default:
    {
        throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::GetParameterDouble] Constitutive law does not have the requested variable");
    }
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::NonlocalDamagePlasticityEngineeringStress::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::BIAXIAL_COMPRESSIVE_STRENGTH:
    {
        this->CheckBiaxialCompressiveStrength(rValue);
        this->mBiaxialCompressiveStrength = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
    {
        this->CheckCompressiveStrength(rValue);
        this->mCompressiveStrength = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::DENSITY:
    {
        this->CheckDensity(rValue);
        this->mRho = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
    {
        this->CheckFractureEnergy(rValue);
        this->mFractureEnergy = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
    {
        this->CheckNonlocalRadius(rValue);
        this->mNonlocalRadius = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
    {
        this->CheckPoissonsRatio(rValue);
        this->mNu = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
    {
        this->CheckTensileStrength(rValue);
        this->mTensileStrength = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
    {
        this->CheckThermalExpansionCoefficient(rValue);
        this->mThermalExpansionCoefficient = rValue;
        this->SetParametersValid();
        break;
    }
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
    {
        this->CheckYoungsModulus(rValue);
        this->mE = rValue;
        this->SetParametersValid();
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::SetParameterDouble] Constitutive law does not have the requested variable");
    }
    }
}


///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::NonlocalDamagePlasticityEngineeringStress::GetType() const
{
    return NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::NonlocalDamagePlasticityEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return true;
    case NuTo::Element::PLANE2D3N:
        return true;
    case NuTo::Element::PLANE2D4N:
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
    case NuTo::Element::ELEMENT1D:
        return true;
    case NuTo::Element::ELEMENT2D:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if density is non negative
//! @param rE ... Young's modulus
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckDensity] The density must not be negative.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check if the nonlocal radius is positive
//! @param rRadius ... nonlocal radius
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckNonlocalRadius(double rRadius) const
{
    if (rRadius <= 0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckNonlocalRadius] Nonlocal radius must be positive.");
    }
}
//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... nonlocal radius
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if compressive strength is positive
//! @param rRadius ... compressive strength
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckCompressiveStrength(double rCompressiveStrength) const
{
    if (rCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckCompressiveStrength] The compressive strength must be a positive value.");
    }
}

//! @brief ... check if biaxial compressive strength is positive
//! @param rBiaxialCompressiveStrength ... biaxial compressive strength
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const
{
    if (rBiaxialCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be a positive value.");
    }
    if (rBiaxialCompressiveStrength <= mCompressiveStrength)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be higher than the uniaxial compressive strength.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}

//! @brief ... check thermal expansion coefficient
//! @param rAlpha ... thermal expansion coefficient
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckThermalExpansionCoefficient(double rAlpha) const
{
}

//! @brief ... calculate the length of the element in plane coordinates (square root of area)
double NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateEquivalentLength2D(const ElementBase* rElement, const Eigen::Matrix<double,4,1>& rStress) const
{
    double l_eq_plane;
    double l_eq_circ;
    double l_element(sqrt(rElement->CalculateArea()));
    l_eq_plane = mNonlocalRadius;
    l_eq_circ  = mNonlocalRadius*rElement->GetSection()->GetThickness()/l_element;

    double factor = 0.5*(rStress[0]-rStress[1]);
    double helpScalar = sqrt(factor*factor+rStress[2]*rStress[2]*0.25);
    double stressPlane;
    double l_eq;
    if ((rStress[0]+rStress[1])*0.5-helpScalar>0)
    {
        // principal plastic strains in plane direction are both in tension -> mixed interpolation
        stressPlane = sqrt(rStress[0]*rStress[0]+rStress[1]*rStress[1]+0.5*rStress[2]*rStress[2]);
    }
    else
    {
        // only the largest component is considered
        stressPlane = (rStress[0]+rStress[1])*0.5+helpScalar;
    }

    if (stressPlane >=rStress[3])
    {
        if (std::abs(rStress[3])>1e-10)
        {
            //plane and zz-direction in tension
            double tan_beta = rStress[3]/stressPlane;
            double x_s = sqrt(1./(1/(l_eq_plane*l_eq_plane)+tan_beta*tan_beta/(l_eq_circ*l_eq_circ)));
            double y_s = x_s * tan_beta;
            l_eq = sqrt(x_s*x_s+y_s*y_s);
        }
        else
        {
            l_eq = l_eq_plane;
        }

    }
    else
    {
        if (std::abs(stressPlane)>1e-10)
        {
            //plane and zz-direction in tension
            double tan_beta = rStress[3]/stressPlane;
            double x_s = sqrt(1./(1/(l_eq_plane*l_eq_plane)+tan_beta*tan_beta/(l_eq_circ*l_eq_circ)));
            double y_s = x_s * tan_beta;
            l_eq = sqrt(x_s*x_s+y_s*y_s);
        }
        else
        {
            l_eq = l_eq_circ;
        }
    }
    return l_eq;
//    return l_eq_plane;
}

//! @brief ... calculate the derivative of the equivalent length of the element in plane coordinates (square root of area) with respect to the plastic strain
double NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateDerivativeEquivalentLength2D(const ElementBase* rElement, const Eigen::Matrix<double,4,1>& rStress,
        Eigen::Matrix<double,4,1>& rdLdStress) const
{
    double l_eq_plane;
    double l_eq_circ;
    double l_element(sqrt(rElement->CalculateArea()));
    l_eq_plane = mNonlocalRadius;
    l_eq_circ  = mNonlocalRadius*rElement->GetSection()->GetThickness()/l_element;
    assert(l_eq_plane>0);
    assert(l_eq_circ>0);

    double factor = 0.5*(rStress[0]-rStress[1]);
    double helpScalar = sqrt(factor*factor+rStress[2]*rStress[2]*0.25);
    double StressPlane;
    double l_eq;
    double d_stress_max_plane_d_stress_xx;
    double d_stress_max_plane_d_stress_yy;
    double d_stress_max_plane_d_stress_xy;

    if ((rStress[0]+rStress[1])*0.5-helpScalar>0)
    {
        // principal plastic strains in plane direction are both in tension -> mixed interpolation
        StressPlane = sqrt(rStress[0]*rStress[0]+rStress[1]*rStress[1]+0.5*rStress[2]*rStress[2]);
        d_stress_max_plane_d_stress_xx = rStress(0)/StressPlane;
        d_stress_max_plane_d_stress_yy = rStress(1)/StressPlane;
        d_stress_max_plane_d_stress_xy = 0.5*rStress(2)/StressPlane;
    }
    else
    {
        // only the largest component is considered
        StressPlane = (rStress[0]+rStress[1])*0.5+helpScalar;

        factor=1./(2.*sqrt((rStress(0)-rStress(1))*(rStress(0)-rStress(1))+rStress(2)*rStress(2)));
        d_stress_max_plane_d_stress_xx = 0.5+(rStress(0)-rStress(1))*factor;
        d_stress_max_plane_d_stress_yy = 0.5+(rStress(1)-rStress(0))*factor;
        d_stress_max_plane_d_stress_xy = rStress(2)*factor;
    }

    if (StressPlane >=rStress[3])
    {
        if (std::abs(rStress[3])>1e-10)
        {
            //plane and zz-direction in tension
            double tan_beta = rStress[3]/StressPlane;
            double x_s = sqrt(1./(1/(l_eq_plane*l_eq_plane)+tan_beta*tan_beta/(l_eq_circ*l_eq_circ)));
            double y_s = x_s * tan_beta;
            l_eq = sqrt(x_s*x_s+y_s*y_s);

            double d_xs_d_tanbeta = -tan_beta*x_s*x_s*x_s/(l_eq_circ*l_eq_circ);
            double d_tan_beta_d_stress_zz = 1./StressPlane;
            double d_tan_beta_d_stress_max_plane = -rStress[3]/(StressPlane*StressPlane);

            Eigen::Matrix<double,4,1> d_xs_dStress;
            d_xs_dStress(0) = d_xs_d_tanbeta *  d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xx;
            d_xs_dStress(1) = d_xs_d_tanbeta *  d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_yy;
            d_xs_dStress(2) = d_xs_d_tanbeta *  d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xy;
            d_xs_dStress(3) = d_xs_d_tanbeta *  d_tan_beta_d_stress_zz;

            Eigen::Matrix<double,4,1> d_ys_dStress;
            d_ys_dStress(0) = tan_beta * d_xs_dStress(0) + d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xx * x_s;
            d_ys_dStress(1) = tan_beta * d_xs_dStress(1) + d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_yy * x_s;
            d_ys_dStress(2) = tan_beta * d_xs_dStress(2) + d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xy * x_s;
            d_ys_dStress(3) = tan_beta * d_xs_dStress(3) + d_tan_beta_d_stress_zz * x_s;

            rdLdStress(0) = x_s/l_eq * d_xs_dStress(0) + y_s/l_eq * d_ys_dStress(0);
            rdLdStress(1) = x_s/l_eq * d_xs_dStress(1) + y_s/l_eq * d_ys_dStress(1);
            rdLdStress(2) = x_s/l_eq * d_xs_dStress(2) + y_s/l_eq * d_ys_dStress(2);
            rdLdStress(3) = x_s/l_eq * d_xs_dStress(3) + y_s/l_eq * d_ys_dStress(3);
            //rLogger<< "[NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateDerivativeEquivalentLength2D rdLdStress1]" << "\n" << rdLdStress << "\n";
        }
        else
        {
            l_eq = l_eq_plane;
            rdLdStress.setZero(4,1);
        }

    }
    else
    {
        if (std::abs(StressPlane)>1e-10)
        {
            //plane and zz-direction in tension
            double tan_beta = rStress[3]/StressPlane;
            double x_s = sqrt(1./(1/(l_eq_plane*l_eq_plane)+tan_beta*tan_beta/(l_eq_circ*l_eq_circ)));
            double y_s = x_s * tan_beta;
            l_eq = sqrt(x_s*x_s+y_s*y_s);

            double d_xs_d_tanbeta = -tan_beta*x_s*x_s*x_s/(l_eq_circ*l_eq_circ);
            double d_tan_beta_d_stress_zz = 1./StressPlane;
            double d_tan_beta_d_stress_max_plane = -rStress[3]/(StressPlane*StressPlane);

            Eigen::Matrix<double,4,1> d_xs_dStress;
            d_xs_dStress(0) = d_xs_d_tanbeta *  d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xx;
            d_xs_dStress(1) = d_xs_d_tanbeta *  d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_yy;
            d_xs_dStress(2) = d_xs_d_tanbeta *  d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xy;
            d_xs_dStress(3) = d_xs_d_tanbeta *  d_tan_beta_d_stress_zz;

            Eigen::Matrix<double,4,1> d_ys_dStress;
            d_ys_dStress(0) = tan_beta * d_xs_dStress(0) + d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xx * x_s;
            d_ys_dStress(1) = tan_beta * d_xs_dStress(1) + d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_yy * x_s;
            d_ys_dStress(2) = tan_beta * d_xs_dStress(2) + d_tan_beta_d_stress_max_plane * d_stress_max_plane_d_stress_xy * x_s;
            d_ys_dStress(3) = tan_beta * d_xs_dStress(3) + d_tan_beta_d_stress_zz;

            rdLdStress(0) = x_s/l_eq * d_xs_dStress(0) + y_s/l_eq * d_ys_dStress(0);
            rdLdStress(1) = x_s/l_eq * d_xs_dStress(1) + y_s/l_eq * d_ys_dStress(1);
            rdLdStress(2) = x_s/l_eq * d_xs_dStress(2) + y_s/l_eq * d_ys_dStress(2);
            rdLdStress(3) = x_s/l_eq * d_xs_dStress(3) + y_s/l_eq * d_ys_dStress(3);
            //rLogger<< "[NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateDerivativeEquivalentLength2D rdLdStress2]" << "\n" << rdLdStress << "\n";
        }
        else
        {
            l_eq = l_eq_circ;
            rdLdStress.setZero(4,1);
        }
    }
    return l_eq;
//    rdLdStress.setZero(4,1);
//    return l_eq_plane;
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
NuTo::Error::eError NuTo::NonlocalDamagePlasticityEngineeringStress::ReturnMapping2D(
        const EngineeringStrain2D& rStrain,
        const double rPrevPlasticStrain[4],
        const EngineeringStrain2D& rPrevTotalStrain,
        Eigen::Matrix<double,4,1>& rStress,
        Eigen::Matrix<double,4,1>& rEpsilonP,
        double& rDeltaEqPlasticStrain,
        Eigen::Matrix<double,4,4>& rdEpsilonPdEpsilon,
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

    //******************************************************************
    //*    F_BETA:    required by DRUCKER-PRAGER yield surface               *
    //******************************************************************
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
    std::vector<Eigen::Matrix<double,4,4> > d2F_d2sigma;
    //! @brief algorithmic modulus * dF_dsigma
    std::vector<Eigen::Matrix<double,4,1> > vectorN;
    //! @brief number of active yield functions
    int numActiveYieldFunctions;

    // for the application of strains in steps, calculate the total strain increment to be applied
    deltaStrain(0) = rStrain[0]-rPrevTotalStrain[0];
    deltaStrain(1) = rStrain[1]-rPrevTotalStrain[1];
    deltaStrain(2) = rStrain[2]-rPrevTotalStrain[2];
    deltaStrain(3) = 0;

    // initialize last plastic strain and last converged stress
    lastPlastStrain << rPrevPlasticStrain[0] , rPrevPlasticStrain[1] ,rPrevPlasticStrain[2] ,rPrevPlasticStrain[3];
    rDeltaEqPlasticStrain = 0;
    //*****************************************************************
    //                  elastic matrix generation                       *
    //*****************************************************************
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

        curTotalStrain(0) = rPrevTotalStrain[0]+cutbackFactorExternal*deltaStrain(0);
        curTotalStrain(1) = rPrevTotalStrain[1]+cutbackFactorExternal*deltaStrain(1);
        curTotalStrain(2) = rPrevTotalStrain[2]+cutbackFactorExternal*deltaStrain(2);
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
                //*************************************************
                //*  thus we have elastic -------------> elastic  *
                //*************************************************
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
            d2F_d2sigma.resize(numYieldSurfaces);

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
                    throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::ReturnMapping2D] programming error - should not happen.");
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
                        if (!YieldSurfaceDruckerPrager2DDerivatives(dF_dsigma[0],&(d2F_d2sigma[0]),trialStress,BETA))
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
                        YieldSurfaceRankine2DRoundedDerivatives(dF_dsigma[1],&(d2F_d2sigma[1]),trialStress);
#ifdef ENABLE_DEBUG
                        rLogger << "dF_dsigma[1] (Rankine)" <<  "\n" << dF_dsigma[1].transpose() << "\n" << "\n";
                        //check second derivative

/*                        double delta=1e-8;
                        rLogger << "d2F_d2sigma ana" << "\n" << d2F_d2sigma[1] << "\n" << "\n";;
                        Eigen::Matrix<double,4,4> dF2_dsigmaTmp;
                        Eigen::Matrix<double,4,4> dF2_dsigmaCD;
                        Eigen::Matrix<double,4,1> dF_dsigmaTmp;
                        for (int count=0; count<4; count++)
                        {
                            trialStress(count,0)+=delta;
                            YieldSurfaceRankine2DRoundedDerivatives(dF_dsigmaTmp,&dF2_dsigmaTmp,trialStress);
                            dF2_dsigmaCD.col(count) = (dF_dsigmaTmp - dF_dsigma[1])*(1./delta);
                            trialStress(count,0)-=delta;
                        }
                        rLogger << "d2F_d2sigma cd" << "\n" << dF2_dsigmaCD << "\n";
                        double max_error((dF2_dsigmaCD-d2F_d2sigma[1]).maxCoeff());
                        rLogger << "max error " << max_error << "\n";

                        if (max_error>1e-6)
                            exit(-1);
*/
#endif
                    }


                    //************************************************************************
                    // residual
                    //************************************************************************
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
                    rLogger << iteration <<" d2F_d2sigma[" << count << "] "<< "\n" << d2F_d2sigma[count] << "\n" << "\n";
#endif
                            hessian+=deltaGamma(count)*d2F_d2sigma[count];
                        }
                    }
#ifdef ENABLE_DEBUG
                    rLogger << iteration <<" hessian" << "\n" << hessian << "\n" << "\n";

                    if (std::abs(hessian.determinant())<toleranceDeterminant)
                    {
                        rLogger << "hessian"<< "\n" << hessian << "\n";
                        rLogger << "trialStress"<< "\n" << trialStress << "\n";
                        rLogger << "yieldConditionFlag " <<  "\n" << yieldConditionFlag.transpose() << "\n" << "\n";
                    }
#endif
                    assert(std::abs(hessian.determinant())>toleranceDeterminant);

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
                            if (std::abs(yieldCondition(count)) > toleranceYieldSurface*f_ct)
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
                            assert(std::abs(matG.determinant())>toleranceDeterminant);
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
                            assert(std::abs(matG.determinant())>toleranceDeterminant);
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
                            //******************************************************************
                            // compute increments for stress
                            //******************************************************************
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


                                //* calculate yield condition and yield condition flag (active or not)
                                //* Drucker Prager
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
            e.AddMessage("[NuTo::NonlocalDamagePlasticityEngineeringStress::ReturnMapping2D] Error performing return mapping procedure.");
            throw;
        }
        catch (...)
        {
               throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::ReturnMapping2D] Error performing return mapping procedure.");
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
        rLogger << "[NuTo::NonlocalDamagePlasticityEngineeringStress::ReturnMapping2D] No convergence can be obtained in the return mapping procedure." << "n";
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
double NuTo::NonlocalDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded(Eigen::Matrix<double,4,1>& rStress, double rFct)const
{
    double value_sum  = 0.5*(rStress(0)+rStress(1));
    double help_scalar = (rStress(0)-rStress(1));
    double value_sqrt = sqrt(help_scalar*help_scalar+4.*rStress(2)*rStress(2));
    double sigma_1(value_sum+0.5*value_sqrt);
    double sigma_2(value_sum-0.5*value_sqrt);
    //* (rounded) Rankine
#ifdef ENABLE_DEBUG
    std::cout << "p1 " << sigma_1 << " p2 " << sigma_2 << " p3 " << rStress[3] << "\n" << "\n";
#endif
    if (rStress[3]<0)
    {
        //sigma_3 is negative
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            if (rStress[3]>sigma_1)
            {
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " all negative f1" << "\n";
#endif
                return rStress[3]-rFct;
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
            return rStress[3]-rFct;
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
                return sqrt(sigma_1*sigma_1 + rStress[3] * rStress[3]) - rFct;
            }
            else
            {
                //sigma_2 is positive
#ifdef ENABLE_DEBUG
                std::cout << "\n" << " f1,f2,f3" << "\n";
#endif
                return sqrt(sigma_1*sigma_1 + sigma_2*sigma_2 + rStress[3] * rStress[3])-rFct;
            }
        }
    }
}

//! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rStress current stress
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
void NuTo::NonlocalDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRoundedDerivatives(Eigen::Matrix<double,4,1>& rdF_dSigma,Eigen::Matrix<double,4,4>* rd2F_d2Sigma,
        Eigen::Matrix<double,4,1>& rStress)const
{
    double value_sum  = 0.5*(rStress(0)+rStress(1));
    double help_diff = (rStress(0)-rStress(1));
    double value_sqrt = sqrt(help_diff*help_diff+4.*rStress(2)*rStress(2));
    double sigma_1(value_sum+0.5*value_sqrt);
    double sigma_2(value_sum-0.5*value_sqrt);

    double factor;
    if (rStress[3]<0)
    {
        //sigma_3 is negative
        if (sigma_1<0)
        {
            // sigma_1 is negative and as a consequence sigma_2 is also negative
            if (rStress[3]>sigma_1)
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
                    throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

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
                    throw MechanicsException("[NuTo::NonlocalDamagePlasticityEngineeringStress::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

                rdF_dSigma(0) = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
                rdF_dSigma(1) = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
                rdF_dSigma(2) = 2*rStress(2)/value_sqrt;
                rdF_dSigma(3) = 0.;

                // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                if (rd2F_d2Sigma!=0)
                {
                    double factor = 1./(value_sqrt*value_sqrt*value_sqrt);
                    (*rd2F_d2Sigma)(0,0) = factor*2.*rStress[2]*rStress[2];
                    (*rd2F_d2Sigma)(0,1) = factor*(-2.)*rStress[2]*rStress[2];
                    (*rd2F_d2Sigma)(0,2) = factor*(2.)*(rStress[1]-rStress[0])*rStress[2];
                    (*rd2F_d2Sigma)(0,3) = 0;
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*2.*rStress[2]*rStress[2];
                    (*rd2F_d2Sigma)(1,2) = factor*(2.)*(rStress[0]-rStress[1])*rStress[2];
                    (*rd2F_d2Sigma)(1,3) = 0;
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = factor*(2.)*(rStress[0]-rStress[1])*(rStress[0]-rStress[1]);
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
                factor = 1./sqrt(rStress[0]*rStress[0]+rStress[1]*rStress[1]+
                                 2*rStress[2]*rStress[2]);
                rdF_dSigma(0) = rStress[0]*factor;
                rdF_dSigma(1) = rStress[1]*factor;
                rdF_dSigma(2) = 2.*rStress[2]*factor;
                rdF_dSigma(3) = 0.;

                if (rd2F_d2Sigma!=0)
                {
                    // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
                    factor = sqrt(rStress[0]*rStress[0]+rStress[1]*rStress[1]+2*rStress[2]*rStress[2]);
                    factor = 1./(factor*factor*factor);

                    (*rd2F_d2Sigma)(0,0) = factor*(rStress[1]*rStress[1]+2.*rStress[2]*rStress[2]);
                    (*rd2F_d2Sigma)(0,1) = -factor*rStress[0]*rStress[1];
                    (*rd2F_d2Sigma)(0,2) = (-2.)*factor*(rStress[0]*rStress[2]);
                    (*rd2F_d2Sigma)(0,3) = 0;
                    (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
                    (*rd2F_d2Sigma)(1,1) = factor*(rStress[0]*rStress[0]+2.*rStress[2]*rStress[2]);
                    (*rd2F_d2Sigma)(1,2) = (-2.)*factor*(rStress[1]*rStress[2]);
                    (*rd2F_d2Sigma)(1,3) = 0;
                    (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
                    (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
                    (*rd2F_d2Sigma)(2,2) = (2.)*factor*(rStress[0]*rStress[0]+rStress[1]*rStress[1]);
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
            throw MechanicsException("[NuTo::Mechanics::NonlocalDamagePlasticityEngineeringStress] Error calculating first derivative of yield function.");
        }

        // std::abs is checked, since of the type of the yield surface changes, the second derivatives are likely to change as well
        if (std::abs(sigma_1)>delta && std::abs(sigma_2)>delta && std::abs(rStress(3))>delta)
        {
            if ((rd2F_d2SigmaCDF-(*rd2F_d2Sigma)).array().abs().maxCoeff()>1e-1 && std::abs(sigma_1)>delta && std::abs(sigma_2)>delta && std::abs(rStress(3))>delta)
            {
                std::cout << "sigmas principal" << sigma_1<< " " << sigma_2 <<  " " << rStress(3) << "\n";
                std::cout << "sigmas " << rStress(0) << " " << rStress(1) <<  " " << rStress(2) <<  " " <<rStress(3) << "\n";
                std::cout << "error second derivatives " << (rd2F_d2SigmaCDF-(*rd2F_d2Sigma)).array().abs().maxCoeff() << "\n";

                std::cout<< "rd2F_d2SigmaCDF " << "\n" << rd2F_d2SigmaCDF << "\n"<< "\n";
                std::cout<< "rd2F_d2Sigma " << "\n" << (*rd2F_d2Sigma) << "\n"<< "\n";
                throw MechanicsException("[NuTo::Mechanics::NonlocalDamagePlasticityEngineeringStress] Error calculating second derivative of yield function.");
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
double NuTo::NonlocalDamagePlasticityEngineeringStress::YieldSurfaceDruckerPrager2D(Eigen::Matrix<double,4,1>& rStress, double rBeta, double rHP)const
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
bool NuTo::NonlocalDamagePlasticityEngineeringStress::YieldSurfaceDruckerPrager2DDerivatives(Eigen::Matrix<double,4,1>& rdF_dSigma,Eigen::Matrix<double,4,4>* rd2F_d2Sigma,
        Eigen::Matrix<double,4,1>& rStress, double rBETA)const
{
    //*******************************************************************
    //*    second invariante                                                 *
    //*******************************************************************
    double invariante_2 = ((rStress(0)-rStress(1))*(rStress(0)-rStress(1))+
                            (rStress(1)-rStress(3))*(rStress(1)-rStress(3))+
                            (rStress(0)-rStress(3))*(rStress(0)-rStress(3)))/6.+
                             rStress(2)*rStress(2);

    if (std::abs(invariante_2)<1e-12)
        return false;
    double factor = 1./(sqrt(invariante_2)*6.);

    rdF_dSigma(0) = factor * (2.*rStress[0]-rStress[1]-rStress[3])+rBETA/3.;
    rdF_dSigma(1) = factor * (2.*rStress[1]-rStress[0]-rStress[3])+rBETA/3.;
    rdF_dSigma(2) = factor * (6.*rStress[2]); /* vector notation from second order tensor */
    rdF_dSigma(3) = factor * (2.*rStress[3]-rStress[0]-rStress[1])+rBETA/3.;

    if (rd2F_d2Sigma!=0)
    {

        //TODO store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
        factor = 1./(invariante_2*sqrt(invariante_2)*12.);
        (*rd2F_d2Sigma)(0,0) = factor * (rStress[1]*rStress[1]-2.*rStress[1]*rStress[3]+rStress[3]*rStress[3]+4.*rStress[2]*rStress[2]);
        (*rd2F_d2Sigma)(0,1) = factor * (-rStress[0]*rStress[1]+rStress[0]*rStress[3]+rStress[1]*rStress[3]-rStress[3]*rStress[3]-2.*rStress[2]*rStress[2]);
        (*rd2F_d2Sigma)(0,2) = factor * (-2.*rStress[2]*(2.*rStress[0]-rStress[1]-rStress[3]));
        (*rd2F_d2Sigma)(0,3) = factor * (rStress[0]*rStress[1]-rStress[0]*rStress[3]-rStress[1]*rStress[1]+rStress[1]*rStress[3]-2.*rStress[2]*rStress[2]);
        (*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
        (*rd2F_d2Sigma)(1,1) = factor * (rStress[0]*rStress[0]-2.*rStress[0]*rStress[3]+rStress[3]*rStress[3]+4.*rStress[2]*rStress[2]);
        (*rd2F_d2Sigma)(1,2) = factor * (2.*rStress[2]*(rStress[0]-2.*rStress[1]+rStress[3]));
        (*rd2F_d2Sigma)(1,3) = factor * (-rStress[0]*rStress[0]+rStress[0]*rStress[1]+rStress[0]*rStress[3]-rStress[1]*rStress[3]-2.*rStress[2]*rStress[2]);
        (*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
        (*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
        (*rd2F_d2Sigma)(2,2) = factor * 4.*(rStress[0]*rStress[0]-rStress[0]*rStress[1]-rStress[0]*rStress[3]+
                                        rStress[1]*rStress[1]-rStress[1]*rStress[3]+rStress[3]*rStress[3]);
        (*rd2F_d2Sigma)(2,3) = factor * 2. *(rStress[0]+rStress[1]-2.*rStress[3])*rStress[2];
        (*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
        (*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
        (*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
        (*rd2F_d2Sigma)(3,3) = factor * (rStress[0]*rStress[0]-2*rStress[0]*rStress[1]+rStress[1]*rStress[1]+4.*rStress[2]*rStress[2]);
    }
    return true;
}


//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
//! @param rLogger stream for the output
void NuTo::NonlocalDamagePlasticityEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
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
void NuTo::NonlocalDamagePlasticityEngineeringStress::CheckParameters()const
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

//! @brief ... calculate the nonlocal equivalente plastic strain of an integration point (scaled with length)
//! @param rElement Element
//! @param rIp integration point
//! @return equivalente plastic strain scaled with length
double NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateNonlocalEquivalentPlasticStrain(const ElementBase* rElement, int rIp, bool& rUnloading)const
{
    double nonlocalKappa(0);

    const std::vector<const NuTo::ElementBase*>& nonlocalElements(rElement->GetNonlocalElements());

    rUnloading = true;
    for (int countNonlocalElement=0; countNonlocalElement<(int)nonlocalElements.size(); countNonlocalElement++)
    {
        const std::vector<double>& weights(rElement->GetNonlocalWeights(rIp,countNonlocalElement));

        assert((int)weights.size()==nonlocalElements[countNonlocalElement]->GetNumIntegrationPoints());

        //Go through all the integration points
        for (int theIP=0; theIP<(int)weights.size(); theIP++)
        {
            const ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain *staticData(nonlocalElements[countNonlocalElement]->GetStaticData(theIP)->AsNonlocalDamagePlasticity2DPlaneStrain());
            nonlocalKappa+=weights[theIP] * staticData->mTmpKappa;
            if (std::abs(staticData->mTmpdEpsilonPdEpsilon[0])>1e-10 ||
                std::abs(staticData->mTmpdEpsilonPdEpsilon[5])>1e-10 ||
                std::abs(staticData->mTmpdEpsilonPdEpsilon[10])>1e-10 ||
                std::abs(staticData->mTmpdEpsilonPdEpsilon[15])>1e-10)
                rUnloading = false;
        }
    }
    //go through all the elements
    return nonlocalKappa;
}

//! @brief ... calculate the nonlocal plastic strain of an integration point
//! @param rElement Element
//! @param rIp integration point
//! @param rNonlocalPlasticStrain nonlocal plastic strain (return value)
//! @return void
void NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateNonlocalPlasticStrain(const ElementBase* rElement, int rIp, double rNonlocalPlasticStrain[4])const
{
    rNonlocalPlasticStrain[0] = 0.;
    rNonlocalPlasticStrain[1] = 0.;
    rNonlocalPlasticStrain[2] = 0.;
    rNonlocalPlasticStrain[3] = 0.;

    const std::vector<const NuTo::ElementBase*>& nonlocalElements(rElement->GetNonlocalElements());

    for (int countNonlocalElement=0; countNonlocalElement<(int)nonlocalElements.size(); countNonlocalElement++)
    {
        const std::vector<double>& weights(rElement->GetNonlocalWeights(rIp,countNonlocalElement));

        assert((int)weights.size()==nonlocalElements[countNonlocalElement]->GetNumIntegrationPoints());

        //Go through all the integration points
        for (int theIP=0; theIP<(int)weights.size(); theIP++)
        {
            const double *tmpPtr(nonlocalElements[countNonlocalElement]->GetStaticData(theIP)->AsNonlocalDamagePlasticity2DPlaneStrain()->mTmpEpsilonP);
            rNonlocalPlasticStrain[0]+=weights[theIP] * tmpPtr[0];
            rNonlocalPlasticStrain[1]+=weights[theIP] * tmpPtr[1];
            rNonlocalPlasticStrain[2]+=weights[theIP] * tmpPtr[2];
            rNonlocalPlasticStrain[3]+=weights[theIP] * tmpPtr[3];
        }
    }
    //go through all the elements
    return;
}


double NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateDamage(double rkappa, double rKappaD)const
{
    double omega = 1.-exp(-rkappa/rKappaD);
    if (omega>MAX_OMEGA)
    {
        omega = MAX_OMEGA;
    }
    return omega;
}

double NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateDerivativeDamage(double rkappa, double rKappaD, double& rdOmegadKappa)const
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


double NuTo::NonlocalDamagePlasticityEngineeringStress::CalculateKappaD()const
{
    return 15.*mFractureEnergy/(16.*mTensileStrength);
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::NonlocalDamagePlasticityEngineeringStress::HaveTmpStaticData() const
{
    return true;
}

//! @brief ... returns true, if a material model has is nonlocal (stiffness is of dynamic size, nonlocal averaging)
//! @return ... see brief explanation
bool NuTo::NonlocalDamagePlasticityEngineeringStress::IsNonlocalModel()const
{
    if (mDamage)
        return true;
    else
        return false;
}
