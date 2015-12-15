// $Id: GradientDamageEngineeringStress.cpp 612 2012-08-13 07:31:23Z unger3 $
// GradientDamageEngineeringStress.cpp
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
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#define MAX_OMEGA 0.999
//#define ENABLE_DEBUG

NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress() :
        ConstitutiveBase(), mRho(0.), mE(0.), mNu(0.), mNonlocalRadius(0.), mNonlocalRadiusParameter(0.), mThermalExpansionCoefficient(0.), mTensileStrength(0.), mCompressiveStrength(0.), mFractureEnergy(
                0.), mDamageLawType(Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING)
{
    SetParametersValid();

#ifdef ENABLE_DEBUG
    std::cout << "NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress debug active" << std::endl;
#else
    std::cout << "NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress debug inactive" << std::endl;
#endif

}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::GradientDamageEngineeringStress::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize GradientDamageEngineeringStress" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
    & BOOST_SERIALIZATION_NVP(mRho)
    & BOOST_SERIALIZATION_NVP(mE)
    & BOOST_SERIALIZATION_NVP(mNu)
    & BOOST_SERIALIZATION_NVP(mNonlocalRadius)
    & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient)
    & BOOST_SERIALIZATION_NVP(mDamageLawType)
    & BOOST_SERIALIZATION_NVP(mDamageLawParameters);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GradientDamageEngineeringStress" << "\n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::GradientDamageEngineeringStress)
#endif // ENABLE_SERIALIZATION

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate1D(ElementBase* rElement, int rIp,
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

    if (section->GetType() != Section::TRUSS)
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] only truss sections are implemented.");

    // calculate engineering strain
    EngineeringStrain1D strain1D;
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] deformation gradient 1d needed to evaluate engineering strain1d.");
    const DeformationGradient1D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D)->second->GetDeformationGradient1D());
    deformationGradient.GetEngineeringStrain(strain1D);

    //Get previous ip_data
    ConstitutiveStaticDataGradientDamage1D *oldStaticData = (rElement->GetStaticData(rIp))->AsGradientDamage1D();

    // nonlocal eq strain
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] nonlocal eq strain needed to evaluate engineering stress 1d.");
    const double nonlocalEqStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)->second->GetNonlocalEqStrain().GetValue(0);

    // last converged kappa
    const double kappaLastConverged = oldStaticData->mKappa;

    // kappa
    const double kappa = std::max(nonlocalEqStrain, kappaLastConverged);

    //calculate damage
    double omega = CalculateDamage(kappa);

    // calculate stress
    EngineeringStress1D stress1D;
    stress1D(0) = mE * strain1D(0);

    // calculate local eq strain and its derivative
    NuTo::FullVector<double, 1> localEqStrainDerivative;
    LocalEqStrain localEqStrain;

    localEqStrain[0] = strain1D[0];
    localEqStrainDerivative[0] = 1.;
    if (strain1D[0] < 0)
    {
        localEqStrain[0] = -strain1D[0];
        localEqStrainDerivative[0] = -1;
    }

    //set this to true, if update is in the map, perform the update after all other outputs have been calculated
    bool performUpdateAtEnd(false);

    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch (itOutput->first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_1D:
        {
            EngineeringStress1D& engineeringStress1D = itOutput->second->GetEngineeringStress1D();
            engineeringStress1D = (1. - omega) * stress1D;
        }
            break;
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
        {
            //this is for the visualize routines
            EngineeringStress3D& engineeringStress3D(itOutput->second->GetEngineeringStress3D());

            engineeringStress3D[0] = (1. - omega) * stress1D[0];
            engineeringStress3D[1] = 0.;
            engineeringStress3D[2] = 0.;
            engineeringStress3D[3] = 0.;
            engineeringStress3D[4] = 0.;
            engineeringStress3D[5] = 0.;
        }
            break;
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
        {
            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
            engineeringStrain3D[0] = strain1D[0];
            engineeringStrain3D[1] = -mNu * strain1D[0];
            engineeringStrain3D[2] = -mNu * strain1D[0];
            engineeringStrain3D[3] = 0.;
            engineeringStrain3D[4] = 0.;
            engineeringStrain3D[5] = 0.;
        }
            break;
        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            LocalEqStrain& localEqStrainOut = itOutput->second->GetLocalEqStrain();
            localEqStrainOut = localEqStrain;
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D:
        {
            ConstitutiveTangentLocal<1, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
            tangent.SetSymmetry(true);
            tangent(0, 0) = (1. - omega) * mE;
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D:
        {
            ConstitutiveTangentLocal<1, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
            tangent.SetSymmetry(true);
            if (nonlocalEqStrain >= kappaLastConverged)
            {
                // loading

                tangent(0, 0) = -CalculateDerivativeDamage(nonlocalEqStrain) * stress1D[0];

            } else
            {
                // unloading
                tangent(0, 0) = 0.;
            }
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D:
        {
            ConstitutiveTangentLocal<1, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
            tangent(0, 0) = localEqStrainDerivative[0];
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_1D:
        {
            ConstitutiveTangentLocal<1, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_1x1());
            tangent.SetSymmetry(true);

            double xi = mNonlocalRadius;
            double dXi = 0;

            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (localEqStrain[0] <= e_xi)
                {
                    double c0 = 0.01;
                    xi = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                    dXi = (mNonlocalRadius - c0) / e_xi;
                }
            }

            double factor = 1. / xi - (localEqStrain[0] - nonlocalEqStrain) / (xi * xi) * dXi;

            tangent(0, 0) = factor * localEqStrainDerivative[0];
        }
            break;
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveTangentLocal<1, 1>& variableNonlocalRadius(itOutput->second->AsConstitutiveTangentLocal_1x1());
            variableNonlocalRadius[0] = mNonlocalRadius;
            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (localEqStrain[0] <= e_xi)
                {
                    double c0 = 0.01;
                    variableNonlocalRadius[0] = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                }
            }
        }
            break;
        case NuTo::Constitutive::Output::DAMAGE:
        {
            itOutput->second->GetDamage().SetDamage(omega);
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate1D] tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            break;
        default:
            throw MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStress::Evaluate1D] output object ") + NuTo::Constitutive::OutputToString(itOutput->first)
                            + std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    //update history variables
    if (performUpdateAtEnd)
    {
        oldStaticData->mKappa = kappa;
    }
    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate2D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    // get section information determining which input on the constitutive level should be used
    //const SectionBase* section(rElement->GetSection());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }

    // calculate engineering strain
    EngineeringStrain2D strain2D;
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate2D] deformation gradient 2d needed to evaluate engineering strain2d.");
    const DeformationGradient2D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D)->second->GetDeformationGradient2D());
    deformationGradient.GetEngineeringStrain(strain2D);

    //Get previous ip_data
    ConstitutiveStaticDataGradientDamage1D *oldStaticData = (rElement->GetStaticData(rIp))->AsGradientDamage1D();

    // nonlocal eq strain
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate2D] nonlocal eq strain needed to evaluate engineering stress 2d.");
    const double nonlocalEqStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)->second->GetNonlocalEqStrain().GetValue(0);

    // last converged kappa
    const double kappaLastConverged = oldStaticData->mKappa;

    // kappa
    const double kappa = std::max(nonlocalEqStrain, kappaLastConverged);

    //calculate damage
    double omega = CalculateDamage(kappa);

//    if (rElement->GetStructure()->ElementGetId(rElement) == 48)
//    {
//        double dKappa = kappa - kappaLastConverged;
////        if (dKappa != 0 and omega > 0.)
//        {
//            for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin();
//                    itOutput != rConstitutiveOutput.end(); itOutput++)
//            {
//                rElement->GetStructure()->GetLogger() << Constitutive::OutputToString(itOutput->first) << "\n";
//            }
//
////            std::cout << "dKappa" << dKappa << std::endl;
//            double dmg0 = CalculateDamage(kappa-dKappa);
//            double dmg1 = omega;
//            double dDmg = CalculateDerivativeDamage(kappa);
//
//            rElement->GetStructure()->GetLogger() << "dDmg " << (dDmg) << "\t kappa " << kappa << "\n";
//
//            rElement->GetStructure()->GetLogger() << "dDamage error [%] " << ((dmg1 - dmg0) / dKappa - dDmg) / dDmg *100 << "\n";
//
//            NuTo::FullVector<double, 3> stress0, stress1, dStress, dStress_CDF;
//            stress0.setZero();
//            stress1.setZero();
//            dStress.setZero();
//            dStress_CDF.setZero();
//            if (nonlocalEqStrain > kappaLastConverged)
//            {
//                double C11, C12, C33;
//                this->CalculateCoefficients3D(C11, C12, C33);
//                // calculate Engineering stress
//                stress0(0) = C11 * strain2D[0] + C12 * strain2D[1];
//                stress0(1) = C11 * strain2D[1] + C12 * strain2D[0];
//                stress0(2) = C33 * strain2D[2] ;
//                stress1 = stress0 * (1-CalculateDamage(kappa-dKappa));
//                stress0 *= (1-omega);
//
//                dStress_CDF = (stress0-stress1)/dKappa;
//
//                dStress(0) = -dDmg * (C11 * strain2D[0] + C12 * strain2D[1]);
//                dStress(1) = -dDmg * (C11 * strain2D[1] + C12 * strain2D[0]);
//                dStress(2) = -dDmg * (C33 * strain2D[2]);
//            }
//
//            rElement->GetStructure()->GetLogger() << dStress << "\n";
//            rElement->GetStructure()->GetLogger() << dStress_CDF << "\n";
//
////            std::cout << "Damage " << omega << std::endl;
////            std::cout << "dDamage " << CalculateDerivativeDamage(kappa) << std::endl;
////            std::cout << "dDamage_CDF " << (dmg1 - dmg0) / dKappa << std::endl;
//        }
//    }

    LocalEqStrain localEqStrain; // = CalculateLocalEqStrain2D(strain2D);
    ConstitutiveTangentLocal<3, 1> localEqStrainTangent; // = CalculateLocalEqStrainTangent2D(strain2D);

    switch (rElement->GetSection()->GetType())
    {
    case Section::PLANE_STRAIN:
        CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStrain(strain2D, localEqStrain, localEqStrainTangent);
        break;
    case Section::PLANE_STRESS:
        CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStress(strain2D, localEqStrain, localEqStrainTangent);
        break;

    default:
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate2D] Invalid type of 2D section behavoir found!!!");
    }

    //set this to true, if update is in the map, perform the update after all other outputs have been calculated
    bool performUpdateAtEnd(false);

    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch (itOutput->first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_2D:
        {
            EngineeringStress2D& engineeringStress2D = itOutput->second->GetEngineeringStress2D();
            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
            {
                double C11, C12, C33;
                this->CalculateCoefficients3D(C11, C12, C33);
                // calculate Engineering stress
                engineeringStress2D[0] = C11 * strain2D[0] + C12 * strain2D[1];
                engineeringStress2D[1] = C11 * strain2D[1] + C12 * strain2D[0];
                engineeringStress2D[2] = C33 * strain2D[2];
                engineeringStress2D *= (1 - omega);
            }
                break;
            case Section::PLANE_STRESS:
                // calculate coefficients of the material matrix
                double C11, C12, C33;
                this->CalculateCoefficients2DPlaneStress(C11, C12, C33);

                // calculate Engineering stress
                engineeringStress2D[0] = C11 * strain2D[0] + C12 * strain2D[1];
                engineeringStress2D[1] = C11 * strain2D[1] + C12 * strain2D[0];
                engineeringStress2D[2] = C33 * strain2D[2];
                engineeringStress2D *= (1 - omega);
                break;
            default:
                throw MechanicsException("[NuTo::GradientDamageEngineeringStress::ENGENEERING_STRESS_2D] Invalid type of 2D section behavoir found!!!");
            }
        }
            break;
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
        {
            //this is for the visualize routines
            EngineeringStress3D& engineeringStress3D(itOutput->second->GetEngineeringStress3D());
            double C11, C12, C33;
            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
                // calculate coefficients of the material matrix
                this->CalculateCoefficients3D(C11, C12, C33);

                // calculate Engineering stress
                engineeringStress3D[0] = C11 * strain2D[0] + C12 * strain2D[1];
                engineeringStress3D[1] = C11 * strain2D[1] + C12 * strain2D[0];
                engineeringStress3D[2] = C12 * (strain2D[0] + strain2D[1]);
                engineeringStress3D[3] = C33 * strain2D[2];
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = 0.;

                engineeringStress3D *= (1 - omega);

                break;
            case Section::PLANE_STRESS:
                // calculate coefficients of the material matrix
                this->CalculateCoefficients2DPlaneStress(C11, C12, C33);

                // calculate Engineering stress
                engineeringStress3D[0] = C11 * strain2D[0] + C12 * strain2D[1];
                engineeringStress3D[1] = C11 * strain2D[1] + C12 * strain2D[0];
                engineeringStress3D[2] = 0.;
                engineeringStress3D[3] = C33 * strain2D[2];
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = 0.;

                engineeringStress3D *= (1 - omega);

                break;
            default:
                throw MechanicsException("[NuTo::GradientDamageEngineeringStress::ENGENEERING_STRESS_3D] Invalid type of 2D section behavoir found!!!");
            }
        }
            break;
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
        {
            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
                engineeringStrain3D[0] = strain2D[0];
                engineeringStrain3D[1] = strain2D[1];
                engineeringStrain3D[2] = 0;
                engineeringStrain3D[3] = strain2D[2];
                engineeringStrain3D[4] = 0.;
                engineeringStrain3D[5] = 0.;
                break;
            case Section::PLANE_STRESS:
                engineeringStrain3D[0] = strain2D[0];
                engineeringStrain3D[1] = strain2D[1];
                engineeringStrain3D[2] = mNu / (mNu - 1.) * (strain2D[0] + strain2D[1]);
                engineeringStrain3D[3] = strain2D[2];
                engineeringStrain3D[4] = 0.;
                engineeringStrain3D[5] = 0.;
                break;
            default:
                throw MechanicsException("[NuTo::GradientDamageEngineeringStress::ENGINEERING_STRAIN_3D] Invalid type of 2D section behavoir found!!!");
            }
        }
            break;
        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            itOutput->second->GetLocalEqStrain() = localEqStrain;
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D:
        {
            ConstitutiveTangentLocal<3, 3>& tangent(itOutput->second->AsConstitutiveTangentLocal_3x3());
            double C11, C12, C33;
            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
                // calculate coefficients of the material matrix

                this->CalculateCoefficients3D(C11, C12, C33);

                // store tangent at the output object
                tangent(0, 0) = C11;
                tangent(1, 0) = C12;
                tangent(2, 0) = 0;

                tangent(0, 1) = C12;
                tangent(1, 1) = C11;
                tangent(2, 1) = 0;

                tangent(0, 2) = 0.;
                tangent(1, 2) = 0.;
                tangent(2, 2) = C33;

                tangent *= (1 - omega);
                break;
            case Section::PLANE_STRESS:
                // calculate coefficients of the material matrix
                this->CalculateCoefficients2DPlaneStress(C11, C12, C33);

                // store tangent at the output object
                tangent(0, 0) = C11;
                tangent(1, 0) = C12;
                tangent(2, 0) = 0;

                tangent(0, 1) = C12;
                tangent(1, 1) = C11;
                tangent(2, 1) = 0;

                tangent(0, 2) = 0.;
                tangent(1, 2) = 0.;
                tangent(2, 2) = C33;

                tangent *= (1 - omega);
                break;
            default:
                throw MechanicsException("[NuTo::GradientDamageEngineeringStress::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D] Invalid type of 2D section behavoir found!!!");
            }

            tangent.SetSymmetry(true);
            break;

        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D:
        {
            ConstitutiveTangentLocal<3, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_3x1());
            tangent.SetSymmetry(false);
            if (nonlocalEqStrain >= kappaLastConverged)
            {
                double damageDerivative = CalculateDerivativeDamage(nonlocalEqStrain);
                // loading
                switch (rElement->GetSection()->GetType())
                {
                case Section::PLANE_STRAIN:
                {
                    double C11, C12, C33;
                    this->CalculateCoefficients3D(C11, C12, C33);
                    // calculate Engineering stress
                    tangent(0) = -damageDerivative * (C11 * strain2D[0] + C12 * strain2D[1]);
                    tangent(1) = -damageDerivative * (C11 * strain2D[1] + C12 * strain2D[0]);
                    tangent(2) = -damageDerivative * (C33 * strain2D[2]);
                }
                    break;
                case Section::PLANE_STRESS:
                    // calculate coefficients of the material matrix
                    double C11, C12, C33;
                    this->CalculateCoefficients2DPlaneStress(C11, C12, C33);

                    // calculate Engineering stress
                    tangent(0) = -damageDerivative * (C11 * strain2D[0] + C12 * strain2D[1]);
                    tangent(1) = -damageDerivative * (C11 * strain2D[1] + C12 * strain2D[0]);
                    tangent(2) = -damageDerivative * (C33 * strain2D[2]);
                    break;
                default:
                    throw MechanicsException("[NuTo::GradientDamageEngineeringStress::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D] Invalid type of 2D section behavoir found!!!");
                }
            } else
            {
                // unloading
                tangent.setZero();
            }
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_2D:
        {
            ConstitutiveTangentLocal<3, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_3x1());
            tangent = localEqStrainTangent;
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_2D:
        {
            ConstitutiveTangentLocal<3, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_3x1());

            double xi = mNonlocalRadius;
            double dXi = 0;

            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (localEqStrain[0] <= e_xi)
                {
                    double c0 = 0.01;
                    xi = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                    dXi = (mNonlocalRadius - c0) / e_xi;
                }
            }

            double factor = 1. / xi - (localEqStrain[0] - nonlocalEqStrain) / (xi * xi) * dXi;

            tangent = factor * localEqStrainTangent;
        }
            break;
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveTangentLocal<1, 1>& variableNonlocalRadius(itOutput->second->AsConstitutiveTangentLocal_1x1());
            variableNonlocalRadius[0] = mNonlocalRadius;
            if (mNonlocalRadiusParameter != 0.)
            {
                double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
                if (localEqStrain[0] <= e_xi)
                {
                    double c0 = 0.01;
                    variableNonlocalRadius[0] = c0 + (mNonlocalRadius - c0) * (localEqStrain[0] / e_xi);
                }
            }
        }
            break;
        case NuTo::Constitutive::Output::DAMAGE:
        {
            itOutput->second->GetDamage().SetDamage(omega);
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate2D] tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            break;
        default:
            throw MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStress::Evaluate2D] output object ") + NuTo::Constitutive::OutputToString(itOutput->first)
                            + std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    //update history variables
    if (performUpdateAtEnd)
    {
        oldStaticData->mKappa = kappa;
    }
    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate3D(ElementBase* rElement, int rIp,
        const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    // get section information determining which input on the constitutive level should be used
    //const SectionBase* section(rElement->GetSection());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }

    // calculate engineering strain
    EngineeringStrain3D strain3D;

    if (rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate2D] deformation gradient 3d needed to evaluate engineering strain3d.");
    const DeformationGradient3D& deformationGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D)->second->GetDeformationGradient3D());
    deformationGradient.GetEngineeringStrain(strain3D);

    //Get previous ip_data
    ConstitutiveStaticDataGradientDamage1D *oldStaticData = (rElement->GetStaticData(rIp))->AsGradientDamage1D();

    // nonlocal eq strain
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN) == rConstitutiveInput.end())
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate3D] nonlocal eq strain needed to evaluate engineering stress 3d.");
    const double nonlocalEqStrain = rConstitutiveInput.find(NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN)->second->GetNonlocalEqStrain().GetValue(0);

    // last converged kappa
    const double kappaLastConverged = oldStaticData->mKappa;

    // kappa
    const double kappa = std::max(nonlocalEqStrain, kappaLastConverged);

    //calculate damage
    double omega = CalculateDamage(kappa);

    LocalEqStrain localEqStrain; // = CalculateLocalEqStrain3D(strain3D);
    ConstitutiveTangentLocal<6, 1> localEqStrainTangent; // = CalculateLocalEqStrainTangent3D(strain3D);

    CalculateLocalEqStrainAndDerivativeModifiedMises3D(strain3D, localEqStrain, localEqStrainTangent);

    //set this to true, if update is in the map, perform the update after all other outputs have been calculated
    bool performUpdateAtEnd(false);

    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch (itOutput->first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_3D:
        {
            EngineeringStress3D& engineeringStress3D = itOutput->second->GetEngineeringStress3D();

            double C11, C12, C44;

            this->CalculateCoefficients3D(C11, C12, C44);
            // calculate Engineering stress
            engineeringStress3D[0] = C11 * strain3D[0] + C12 * strain3D[1] + C12 * strain3D[2];
            engineeringStress3D[1] = C11 * strain3D[1] + C12 * strain3D[0] + C12 * strain3D[2];
            engineeringStress3D[2] = C11 * strain3D[2] + C12 * strain3D[0] + C12 * strain3D[1];
            engineeringStress3D[3] = C44 * strain3D[3];
            engineeringStress3D[4] = C44 * strain3D[4];
            engineeringStress3D[5] = C44 * strain3D[5];

            engineeringStress3D *= (1 - omega);
        }
            break;
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D:
        {
            EngineeringStrain3D& engineeringStrain3D(itOutput->second->GetEngineeringStrain3D());
            engineeringStrain3D[0] = strain3D[0];
            engineeringStrain3D[1] = strain3D[1];
            engineeringStrain3D[2] = strain3D[2];
            engineeringStrain3D[3] = strain3D[3];
            engineeringStrain3D[4] = strain3D[4];
            engineeringStrain3D[5] = strain3D[5];
        }
            break;
        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            itOutput->second->GetLocalEqStrain() = localEqStrain;
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D:
        {
            ConstitutiveTangentLocal<6, 6>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x6());
            double C11, C12, C44;
            // calculate coefficients of the material matrix

            this->CalculateCoefficients3D(C11, C12, C44);

            // store tangent at the output object
            tangent(0, 0) = C11;
            tangent(1, 0) = C12;
            tangent(2, 0) = C12;
            tangent(3, 0) = 0;
            tangent(4, 0) = 0;
            tangent(5, 0) = 0;

            tangent(0, 1) = C12;
            tangent(1, 1) = C11;
            tangent(2, 1) = C12;
            tangent(3, 1) = 0;
            tangent(4, 1) = 0;
            tangent(5, 1) = 0;

            tangent(0, 2) = C12;
            tangent(1, 2) = C12;
            tangent(2, 2) = C11;
            tangent(3, 2) = 0;
            tangent(4, 2) = 0;
            tangent(5, 2) = 0;

            tangent(0, 3) = 0;
            tangent(1, 3) = 0;
            tangent(2, 3) = 0;
            tangent(3, 3) = C44;
            tangent(4, 3) = 0;
            tangent(5, 3) = 0;

            tangent(0, 4) = 0;
            tangent(1, 4) = 0;
            tangent(2, 4) = 0;
            tangent(3, 4) = 0;
            tangent(4, 4) = C44;
            tangent(5, 4) = 0;

            tangent(0, 5) = 0;
            tangent(1, 5) = 0;
            tangent(2, 5) = 0;
            tangent(3, 5) = 0;
            tangent(4, 5) = 0;
            tangent(5, 5) = C44;

            tangent *= (1 - omega);
            tangent.SetSymmetry(true);
        }
            break;
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_3D:
        {
            ConstitutiveTangentLocal<6, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x1());
            if (nonlocalEqStrain >= kappaLastConverged)
            {
                double damageDerivative = CalculateDerivativeDamage(nonlocalEqStrain);
                // loading
                double C11, C12, C44;
                this->CalculateCoefficients3D(C11, C12, C44);
                // calculate Engineering stress
                tangent(0) = -damageDerivative * (C11 * strain3D[0] + C12 * strain3D[1] + C12 * strain3D[2] );
                tangent(1) = -damageDerivative * (C11 * strain3D[1] + C12 * strain3D[0] + C12 * strain3D[2] );
                tangent(2) = -damageDerivative * (C11 * strain3D[2] + C12 * strain3D[0] + C12 * strain3D[1] );

                tangent(3) = -damageDerivative * (C44 * strain3D[3]);
                tangent(4) = -damageDerivative * (C44 * strain3D[4]);
                tangent(5) = -damageDerivative * (C44 * strain3D[5]);

            } else
            {
                // unloading
                tangent.setZero();
            }
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_3D:
        {
            ConstitutiveTangentLocal<6, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x1());
            tangent = localEqStrainTangent;
        }
            break;
        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D:
        {
            ConstitutiveTangentLocal<6, 1>& tangent(itOutput->second->AsConstitutiveTangentLocal_6x1());

            double xi = mNonlocalRadius;
            double dXi = 0;

            if (mNonlocalRadiusParameter != 0.)
            {
                throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate3D] D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D not implemented for 3D for mNonlocalRadiusParameter != 0.");
            }

            double factor = 1. / xi - (localEqStrain[0] - nonlocalEqStrain) / (xi * xi) * dXi;

            tangent = factor * localEqStrainTangent;

        }
            break;
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveTangentLocal<1, 1>& variableNonlocalRadius(itOutput->second->AsConstitutiveTangentLocal_1x1());
            variableNonlocalRadius[0] = mNonlocalRadius;
            if (mNonlocalRadiusParameter != 0.)
            {
                throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate3D] NONLOCAL_PARAMETER_XI not implemented for 3D for mNonlocalRadiusParameter != 0.");
            }
        }
            break;
        case NuTo::Constitutive::Output::DAMAGE:
        {
            itOutput->second->GetDamage().SetDamage(omega);
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException("[NuTo::GradientDamageEngineeringStress::Evaluate3D] tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            break;
        default:
            throw MechanicsException(
                    std::string("[NuTo::GradientDamageEngineeringStress::Evaluate3D] output object ") + NuTo::Constitutive::OutputToString(itOutput->first)
                            + std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    //update history variables
    if (performUpdateAtEnd)
    {
        oldStaticData->mKappa = kappa;
    }
    return Error::SUCCESSFUL;
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage1D();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage1D();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage1D();
}

// calculate coefficients of the material matrix
void NuTo::GradientDamageEngineeringStress::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE / ((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE / (2 * (1.0 + this->mNu));
}

// calculate coefficients of the material matrix
void NuTo::GradientDamageEngineeringStress::CalculateCoefficients2DPlaneStress(double& C11, double& C12, double& C33) const
{
    double factor = this->mE / (1.0 - (this->mNu * this->mNu));
    C11 = factor;
    C12 = factor * this->mNu;
    C33 = factor * 0.5 * (1.0 - this->mNu);
}
// parameters /////////////////////////////////////////////////////////////

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::GradientDamageEngineeringStress::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        return mCompressiveStrength;
    case Constitutive::eConstitutiveParameter::DENSITY:
        return this->mRho;
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        return mFractureEnergy;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        return this->mNonlocalRadius;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
        return mNonlocalRadiusParameter;
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
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::GetParameterDouble] Constitutive law does not have the requested variable");
    }
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::GradientDamageEngineeringStress::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
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
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
    {
        this->CheckNonlocalRadiusParameter(rValue);
        mNonlocalRadiusParameter = rValue;
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
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::SetParameterDouble] Constitutive law does not have the requested variable");
    }
    }
}


//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::GradientDamageEngineeringStress::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
    {
        NuTo::FullVector<double, Eigen::Dynamic> damageLaw(mDamageLawParameters.rows() + 1);

        damageLaw[0] = static_cast<double>(mDamageLawType);
        damageLaw.SetBlock(1, 0, mDamageLawParameters);
        return damageLaw;
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::GetParameterFullVectorDouble] Constitutive law does not have the requested variable");
    }
    }
}


//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::GradientDamageEngineeringStress::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
    {
        this->CheckDamageLaw(rValue);
        mDamageLawType = static_cast<int>(rValue[0]);
        int numDamageLawParameters = rValue.rows() - 1;
        if (numDamageLawParameters > 0)
            mDamageLawParameters = rValue.GetBlock(1, 0, numDamageLawParameters, 1);

        this->SetParametersValid();
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::GradientDamageEngineeringStress::SetParameterFullVectorDouble] Constitutive law does not have the requested variable");
    }
    }
}


///////////////////////////////////////////////////////////////////////////

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::GradientDamageEngineeringStress::GetType() const
{
    return NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::GradientDamageEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::ELEMENT1D:
        return true;
    case NuTo::Element::ELEMENT2D:
        return true;
    case NuTo::Element::ELEMENT3D:
//        return true;
    case NuTo::Element::BOUNDARYELEMENT1D:
        return true;
    case NuTo::Element::BOUNDARYELEMENT2D:
        return true;
//    case NuTo::Element::BOUNDARYELEMENT3D:
//        return true;
    default:
        return false;
    }
}

//! @brief ... check if density is non negative
//! @param rE ... Young's modulus
void NuTo::GradientDamageEngineeringStress::CheckDensity(double rRho) const
{
    if (rRho < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckDensity] The density must not be negative.");
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::GradientDamageEngineeringStress::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::GradientDamageEngineeringStress::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check if the nonlocal radius is positive
//! @param rRadius ... nonlocal radius
void NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius(double rRadius) const
{
    if (rRadius <= 0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius] Nonlocal radius must be positive.");
    }
}

//! @brief ... check if the nonlocal radius parameter is greater than 1
//! @param rRadius ... nonlocal radius
void NuTo::GradientDamageEngineeringStress::CheckNonlocalRadiusParameter(double rRadiusParameter) const
{
    if (rRadiusParameter <= 1 and rRadiusParameter != 0.)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckNonlocalRadius] Nonlocal radius parameter must be greater than 1. or 0.");
    }
}

//! @brief ... check thermal expansion coefficient
//! @param rAlpha ... thermal expansion coefficient
void NuTo::GradientDamageEngineeringStress::CheckThermalExpansionCoefficient(double rAlpha) const
{
}

//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... nonlocal radius
void NuTo::GradientDamageEngineeringStress::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if compressive strength is positive
//! @param rRadius ... compressive strength
void NuTo::GradientDamageEngineeringStress::CheckCompressiveStrength(double rCompressiveStrength) const
{
    if (rCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckCompressiveStrength] The compressive strength must be a positive value.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::GradientDamageEngineeringStress::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}

//! @brief ... check damage law parameters
//! @param rDamageLawParameters ... damage law parameters
void NuTo::GradientDamageEngineeringStress::CheckDamageLaw(const NuTo::FullVector<double, Eigen::Dynamic>& rDamageLaw) const
{
    int damageLawType = static_cast<int>(rDamageLaw[0]);
    int numDamageLawParameters = rDamageLaw.rows() - 1;

    switch (damageLawType)
    {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
        break;
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        if (numDamageLawParameters != 3)
            throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStress::CheckDamageLaw] ") + std::string("ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH needs exactly 3 parameters."));

    }
        break;
    default:
    {
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStress::CheckDamageLaw] ") + std::string("The required damage law is not implemented. "));
        break;
    }

    }

}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
//! @param rLogger stream for the output
void NuTo::GradientDamageEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus      : " << this->mE << "\n";
    rLogger << "    Poisson's ratio      : " << this->mNu << "\n";
    rLogger << "    nonlocal radius      : " << this->mNonlocalRadius << "\n";
    rLogger << "    damage law params    : " << this->mDamageLawParameters << "\n";
    rLogger << "    thermal expansion coeff : " << this->mThermalExpansionCoefficient << "\n";
}

// check parameters
void NuTo::GradientDamageEngineeringStress::CheckParameters() const
{
    this->CheckDensity(mRho);
    this->CheckYoungsModulus(mE);
    this->CheckPoissonsRatio(mNu);
    this->CheckNonlocalRadius(mNonlocalRadius);
    this->CheckTensileStrength(mTensileStrength);
    this->CheckCompressiveStrength(mCompressiveStrength);
    this->CheckFractureEnergy(mFractureEnergy);
    this->CheckDamageLaw(GetParameterFullVectorDouble(Constitutive::eConstitutiveParameter::DAMAGE_LAW));
    this->CheckThermalExpansionCoefficient(mThermalExpansionCoefficient);
}

NuTo::LocalEqStrain NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrain2D(const NuTo::EngineeringStrain2D& rStrain2D) const
{
    NuTo::LocalEqStrain localEqStrain;
    // calculate principal strains e1 and e2

    /**                          _____________________
     /    2              2
     exx   eyy      /  exy    /exx   eyy\
   e1 =   --- + --- +   /   ---- + |--- - ---|
     2     2    \/     4     \ 2     2 /

     A     +               B
     **/
    double A = (rStrain2D(0) + rStrain2D(1)) / 2.;
    double B = std::sqrt(std::pow(rStrain2D(0) - rStrain2D(1), 2.) / 4. + std::pow(rStrain2D(2), 2.) / 4.);

    // macaulay brackets of principal strains
    double e1 = std::max(A + B, 0.);
    double e2 = std::max(A - B, 0.);

    localEqStrain(0) = std::sqrt(e1 * e1 + e2 * e2);
    return localEqStrain;
}

NuTo::LocalEqStrain NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrain3D(const NuTo::EngineeringStrain3D& rStrain3D) const
{
    throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrain3D] not implemented");
}

NuTo::ConstitutiveTangentLocal<3, 1> NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent2D(const NuTo::EngineeringStrain2D& rStrain2D) const
{
    NuTo::ConstitutiveTangentLocal<3, 1> tangent;
    double A = (rStrain2D(0) + rStrain2D(1)) / 2.;
    double B = 0.5 * std::sqrt(std::pow(rStrain2D(0) - rStrain2D(1), 2.) + std::pow(rStrain2D(2), 2.));

    // macaulay brackets of principal strains
    double e1 = std::max(A + B, 0.);
    double e2 = std::max(A - B, 0.);

    if (e1 + e2 == 0)
    {
        tangent.setZero();
        return tangent;
    }
    if (e1 != 0 and e2 != 0)
    {
        double eMazar = std::sqrt(e1 * e1 + e2 * e2);
        tangent(0) = rStrain2D(0) / eMazar;
        tangent(1) = rStrain2D(1) / eMazar;
        tangent(2) = rStrain2D(2) / eMazar * 0.5;
        return tangent;
    }
    if (e1 != 0)
    {
        tangent(0) = (e1 - rStrain2D[1]) / (2 * B);
        tangent(1) = (e1 - rStrain2D[0]) / (2 * B);
        tangent(2) = 0.5 * rStrain2D[2] / (2 * B);
        return tangent;
    }
    if (e2 != 0)
    {
        tangent(0) = (rStrain2D[0] - e2) / (2 * B);
        tangent(1) = (rStrain2D[1] - e2) / (2 * B);
        tangent(2) = -0.5 * rStrain2D[2] / (2 * B);
        return tangent;
    }

//    std::cout << e1 << std::endl;
//    std::cout << e2 << std::endl;

    throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent2D] error");
}

NuTo::ConstitutiveTangentLocal<6, 1> NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent3D(const NuTo::EngineeringStrain3D& rStrain3D) const
{
    throw NuTo::MechanicsException("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent3D] not implemented");
}

//! @brief calculates the local eq strain and its derivatives for the modified mises model
//! @param rStrain2D ... 2d strain
//! @param rLocalEqStrain ... local eq strain
//! @param rLocalEqStrainTangent ... d local eq strain / d strain 2D
void NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStrain(const NuTo::EngineeringStrain2D& rStrain2D, NuTo::LocalEqStrain& rLocalEqStrain,
        NuTo::ConstitutiveTangentLocal<3, 1>& rLocalEqStrainTangent) const
{
    double k = mCompressiveStrength / mTensileStrength;

    double K1 = (k - 1.) / (2. * k * (1 - 2 * mNu));
    double K2 = 3 / (k * (1 + mNu) * (1 + mNu));

    double I1 = rStrain2D[0] + rStrain2D[1];
    double J2 = 1. / 3. * (rStrain2D[0] * rStrain2D[0] + rStrain2D[1] * rStrain2D[1] - rStrain2D[0] * rStrain2D[1]) + 0.25 * rStrain2D[2] * rStrain2D[2];

    double A = K1 * K1 * I1 * I1 + K2 * J2;
    if (A == 0) // prevent division by 0
    {
        rLocalEqStrain[0] = 0;

        rLocalEqStrainTangent[0] = K1;
        rLocalEqStrainTangent[1] = K1;
        rLocalEqStrainTangent[2] = 0;
        return;
    }

    rLocalEqStrain[0] = K1 * I1 + std::sqrt(A);

    double dJ2dexx = 1. / 3. * (2 * rStrain2D[0] - rStrain2D[1]);
    double dJ2deyy = 1. / 3. * (2 * rStrain2D[1] - rStrain2D[0]);
    double dJ2dgxy = 0.5 * rStrain2D[2];

    rLocalEqStrainTangent[0] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2dexx);
    rLocalEqStrainTangent[1] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2deyy);
    rLocalEqStrainTangent[2] = 1. / (2 * std::sqrt(A)) * (K2 * dJ2dgxy);

}

//! @brief calculates the local eq strain and its derivatives for the modified mises model
//! @param rStrain3D ... 3d strain
//! @param rLocalEqStrain ... local eq strain
//! @param rLocalEqStrainTangent ... d local eq strain / d strain 3D
void NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainAndDerivativeModifiedMises3D(const NuTo::EngineeringStrain3D& rStrain3D, NuTo::LocalEqStrain& rLocalEqStrain,
        NuTo::ConstitutiveTangentLocal<6, 1>& rLocalEqStrainTangent) const
{
    const double k = mCompressiveStrength / mTensileStrength;
    const double K1 = (k - 1.) / (2. * k * (1 - 2 * mNu));
    const double K2 = 3. / (k * (1 + mNu) * (1 + mNu));

    const double eps_xx = rStrain3D[0];
    const double eps_yy = rStrain3D[1];
    const double eps_zz = rStrain3D[2];
    const double eps_xy = 0.5*rStrain3D[3];
    const double eps_yz = 0.5*rStrain3D[4];
    const double eps_zx = 0.5*rStrain3D[5];

    double I1 = eps_xx + eps_yy + eps_zz;
    double J2 = 1. / 6. * (std::pow(eps_xx - eps_yy, 2) + std::pow(eps_yy - eps_zz, 2) + std::pow(eps_zz - eps_xx, 2))
            + std::pow(eps_xy, 2)+ std::pow(eps_yz, 2)+ std::pow(eps_zx, 2);

    double A = K1 * K1 * I1 * I1 + K2 * J2;

    if (A == 0) // prevent division by 0
    {
        rLocalEqStrain[0] = 0;

        rLocalEqStrainTangent[0] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[1] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[2] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[3] = 0;
        rLocalEqStrainTangent[4] = 0;
        rLocalEqStrainTangent[5] = 0;

        return;
    }

    rLocalEqStrain[0] = K1 * I1 + std::sqrt(A);

    const double dJ2dexx = 1. / 3. * (2*eps_xx - eps_yy - eps_zz);
    const double dJ2deyy = 1. / 3. * (2*eps_yy - eps_xx - eps_zz);
    const double dJ2dezz = 1. / 3. * (2*eps_zz - eps_xx - eps_yy);
    const double dJ2dgxy = eps_xy;
    const double dJ2dgyz = eps_yz;
    const double dJ2dgzx = eps_zx;

    rLocalEqStrainTangent[0] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2dexx);
    rLocalEqStrainTangent[1] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2deyy);
    rLocalEqStrainTangent[2] = K1 + 1. / (2 * std::sqrt(A)) * (2 * K1 * K1 * I1 + K2 * dJ2dezz);
    rLocalEqStrainTangent[3] = 1. / (2 * std::sqrt(A)) * K2 * dJ2dgxy;
    rLocalEqStrainTangent[4] = 1. / (2 * std::sqrt(A)) * K2 * dJ2dgyz;
    rLocalEqStrainTangent[5] = 1. / (2 * std::sqrt(A)) * K2 * dJ2dgzx;

}

//! @brief calculates the local eq strain and its derivatives for the modified mises model in plane stress state
//! @param rStrain2D ... 2d strain
//! @param rLocalEqStrain ... local eq strain
//! @param rLocalEqStrainTangent ... d local eq strain / d strain 2D
void NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStress(const NuTo::EngineeringStrain2D& rStrain2D, NuTo::LocalEqStrain& rLocalEqStrain,
        NuTo::ConstitutiveTangentLocal<3, 1>& rLocalEqStrainTangent) const
{
    const double k = mCompressiveStrength / mTensileStrength;
    const double K1 = (k - 1.) / (2. * k * (1 - 2 * mNu));
    const double K2 = 3. / (k * (1 + mNu) * (1 + mNu));

    const double eps_xx = rStrain2D[0];
    const double eps_yy = rStrain2D[1];
    const double eps_zz = mNu / (mNu - 1) * (rStrain2D[0] + rStrain2D[1]);
    const double eps_xy = 0.5 * rStrain2D[2];

    double I1 = eps_xx + eps_yy + eps_zz;
    double J2 = 1. / 6. * (std::pow(eps_xx - eps_yy, 2) + std::pow(eps_yy - eps_zz, 2) + std::pow(eps_zz - eps_xx, 2)) + eps_xy * eps_xy;

    double A = K1 * K1 * I1 * I1 + K2 * J2;

    if (A == 0) // prevent division by 0
    {
        rLocalEqStrain[0] = 0;

        rLocalEqStrainTangent[0] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[1] = (1 + mNu / (mNu - 1)) * K1;
        rLocalEqStrainTangent[2] = 0;

        return;
    }

    rLocalEqStrain[0] = K1 * I1 + std::sqrt(A);

    double dJ2dexx = 1. / 3. * (2. * eps_xx - eps_yy - 2. * eps_zz + 2 * mNu / (mNu - 1) * eps_zz);
    double dJ2deyy = 1. / 3. * (2. * eps_yy - eps_xx - 2. * eps_zz + 2 * mNu / (mNu - 1) * eps_zz);
    double dJ2dgxy = eps_xy;

    rLocalEqStrainTangent[0] = (1 + mNu / (mNu - 1)) * K1 + 1. / (2. * std::sqrt(A)) * (2. * (1 + mNu / (mNu - 1)) * K1 * K1 * I1 + K2 * dJ2dexx);
    rLocalEqStrainTangent[1] = (1 + mNu / (mNu - 1)) * K1 + 1. / (2. * std::sqrt(A)) * (2. * (1 + mNu / (mNu - 1)) * K1 * K1 * I1 + K2 * dJ2deyy);
    rLocalEqStrainTangent[2] = 1. / (2. * std::sqrt(A)) * (K2 * dJ2dgxy);

}
double NuTo::GradientDamageEngineeringStress::CalculateDamage(double rKappa) const
{
    double omega = 0;

    double e_0 = mTensileStrength / mE;
    double e_f = mFractureEnergy / mTensileStrength;
    double e_c = 2 * e_f; // or something

    switch (mDamageLawType)
    {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
    {
        if (rKappa > e_0)
        {
            omega = 1 - e_0 / rKappa;
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
    {
        if (rKappa > e_0)
        {
            omega = e_c / rKappa * (rKappa - e_0) / (e_c - e_0);
            omega = std::min(omega, MAX_OMEGA);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
    {
        if (rKappa > e_0)
        {
            omega = 1 - e_0 / rKappa * exp((e_0 - rKappa) / e_f);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
    {
        if (rKappa > e_0)
        {
            omega = 1 - e_0 / rKappa * (1 - MAX_OMEGA + MAX_OMEGA * exp((e_0 - rKappa) / e_f));
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        double alpha = mDamageLawParameters[0];
        double beta = mDamageLawParameters[1];
        double gamma = mDamageLawParameters[2];

        double e_Dev_e0 = rKappa / e_0;
        double e0_Dev_e = e_0 / rKappa;

        if (e_Dev_e0 <= alpha)
            break;

        if (e_Dev_e0 < beta)
        {
            // pre peak polynomial smoothing
            omega = 1. - e0_Dev_e * ((1 - alpha) * pow((e_Dev_e0 - beta) / (beta - alpha), 3.) + 1);
            break;
        }

        double term = 2. * e_f / (e_0 * (gamma - beta));
        double A = -1. / (1 + term);
        if (e_Dev_e0 < gamma)
        {
            // post peak polynomial smoothing
            omega = 1. - e0_Dev_e * (pow((e_Dev_e0 - beta) / (gamma - beta), 2.) * A + 1);
            break;
        }
        double e_s = e_0 * gamma + e_f * log(-A * term);

        // else
        omega = 1. - e0_Dev_e * exp((e_s - rKappa) / e_f);

        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
    {
        if (rKappa > e_0)
        {
            if (rKappa < e_c)
            {
                double kappa_scaled = (rKappa - e_0) / (e_c - e_0);
                omega = 1 - e_0 / rKappa * (2 * kappa_scaled * kappa_scaled * kappa_scaled - 3 * kappa_scaled * kappa_scaled + 1);
            } else
            {
                omega = 1.;
            }
        }
        break;
    }
    default:
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStress::CalculateDamage] ") + std::string("The required damage law is not implemented. "));

        break;
    }

    return omega;
}

double NuTo::GradientDamageEngineeringStress::CalculateDerivativeDamage(double rKappa) const
{
    double DomegaDkappa = 0;

    double e_0 = mTensileStrength / mE;
    double e_f = mFractureEnergy / mTensileStrength;
    double e_c = 2 * e_f; // or something

    switch (mDamageLawType)
    {
    case Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = e_0 / (rKappa * rKappa);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING:
    {
        double termA = e_c / MAX_OMEGA / (e_c - e_0);
        double kappa_max = termA * e_0 / (termA - 1);

        if (rKappa > kappa_max)
            std::cout << CalculateDamage(kappa_max) << std::endl;

        if (rKappa > e_0 && rKappa < kappa_max)
        {
            DomegaDkappa = e_c * e_0 / (rKappa * rKappa * (e_c - e_0));
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = e_0 / rKappa * (1 / rKappa + 1 / e_f) * exp((e_0 - rKappa) / e_f);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = e_0 / rKappa * ((1 / rKappa + 1 / e_f) * MAX_OMEGA * exp((e_0 - rKappa) / e_f) + (1 - MAX_OMEGA) / rKappa);
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH:
    {
        double alpha = mDamageLawParameters[0];
        double beta = mDamageLawParameters[1];
        double gamma = mDamageLawParameters[2];

        double e_Dev_e0 = rKappa / e_0;
        double e0_dev_e2 = e_0 / pow(rKappa, 2.);

        if (e_Dev_e0 <= alpha)
            break;

        if (e_Dev_e0 < beta)
        {
            // pre peak polynomial smoothing
            DomegaDkappa = e0_dev_e2 * (1. + (1. - alpha) / pow(beta - alpha, 3.) * pow(e_Dev_e0 - beta, 2.) * (-2. * e_Dev_e0 - beta));
            break;
        }

        double term = 2. * e_f / (e_0 * (gamma - beta));
        double A = -1. / (1 + term);
        if (e_Dev_e0 < gamma)
        {
            // post peak polynomial smoothing
            DomegaDkappa = e0_dev_e2 * (1. + (beta * beta - e_Dev_e0 * e_Dev_e0) / ((gamma - beta) * (gamma - beta - 2. * e_f / e_0)));
            break;
        }
        double e_s = e_0 * gamma + e_f * log(-A * term);

        // else
        DomegaDkappa = e0_dev_e2 * (rKappa / e_f + 1.) * exp((e_s - rKappa) / e_f);

        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE:
    {
        if (rKappa > e_0 && rKappa < e_c)
        {
            double kappa_scaled = (rKappa - e_0) / (e_c - e_0);
            DomegaDkappa = -6 * e_0 / rKappa / (e_c - e_0) * (kappa_scaled * kappa_scaled - kappa_scaled)
                    + e_0 / (rKappa * rKappa) * (2 * kappa_scaled * kappa_scaled * kappa_scaled - 3 * kappa_scaled * kappa_scaled + 1);

        }
        break;
    }
    default:
        throw NuTo::MechanicsException(std::string("[NuTo::GradientDamageEngineeringStress::CalculateDerivativeDamage] ") + std::string("The required damage law is not implemented. "));

        break;
    }

    return DomegaDkappa;
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::GradientDamageEngineeringStress::HaveTmpStaticData() const
{
    return false;
}
