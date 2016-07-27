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

#include "nuto/mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "nuto/mechanics/constitutive/inputoutput/EquivalentStrain.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataGradientDamage.h"

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/elements/IpDataStaticDataBase.h"

#define MAX_OMEGA 0.999
//#define ENABLE_DEBUG

NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress() :
        ConstitutiveBase(),
        mRho(0.),
        mE(0.),
        mNu(0.),
        mNonlocalRadius(0.),
        mNonlocalRadiusParameter(0.),
        mThermalExpansionCoefficient(0.),
        mTensileStrength(0.),
        mCompressiveStrength(0.),
        mFractureEnergy(0.),
        mDamageLawType(Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING)
{}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
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
       & BOOST_SERIALIZATION_NVP(mNonlocalRadiusParameter)
       & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient)
       & BOOST_SERIALIZATION_NVP(mTensileStrength)
       & BOOST_SERIALIZATION_NVP(mCompressiveStrength)
       & BOOST_SERIALIZATION_NVP(mFractureEnergy)
       & BOOST_SERIALIZATION_NVP(mDamageLawType);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GradientDamageEngineeringStress \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::GradientDamageEngineeringStress)
#endif // ENABLE_SERIALIZATION

NuTo::ConstitutiveInputMap NuTo::GradientDamageEngineeringStress::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    // Always allocate these two guys. For a dof separated approach, this could be handled more efficiently.
    constitutiveInputMap[Constitutive::Input::ENGINEERING_STRAIN];
    constitutiveInputMap[Constitutive::Input::NONLOCAL_EQ_STRAIN];

    return constitutiveInputMap;
}

NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate1D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    if (rElement->GetSection()->GetType() != Section::TRUSS)
        throw MechanicsException(__PRETTY_FUNCTION__, "only truss sections are implemented.");

    // get constitutive inputs
    const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain1D();
    const auto& nonlocalEqStrain = *rConstitutiveInput.at(Constitutive::Input::NONLOCAL_EQ_STRAIN);

    auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<1>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

    ConstitutiveStaticDataGradientDamage currentStaticData = GetCurrentStaticData(*rElement, rIp, rConstitutiveInput);
    double omega = CalculateDamage(currentStaticData.GetKappa());

    EquivalentStrainModifiedMises<1> eeq(elasticEngineeringStrain, mCompressiveStrength/mTensileStrength, mNu);
    double localEqStrain = eeq.Get();

    bool performUpdateAtEnd = false;

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            engineeringStress[0] = (1. - omega) * mE *  elasticEngineeringStrain[0];
            break;
        }

        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& localEqStrainOut = *itOutput.second;
            localEqStrainOut.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            localEqStrainOut[0] = localEqStrain;
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1,1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = (1. - omega) * mE;
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);
            if (nonlocalEqStrain[0] == currentStaticData.GetKappa())
            {
                // loading
                tangent(0, 0) = -CalculateDerivativeDamage(nonlocalEqStrain[0]) * mE * elasticEngineeringStrain[0];

            } else
            {
                // unloading
                tangent(0, 0) = 0.;
            }
            break;
        }

        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);
            tangent[0] = eeq.GetDerivative()[0];
            break;
        }

        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent[0] = CalculateLocalEqStrainXiFactor(localEqStrain, nonlocalEqStrain[0]) * eeq.GetDerivative()[0];
            break;
        }
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveIOBase& xi = *itOutput.second;
            xi.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            xi[0] = CalculateXi(localEqStrain);
            break;
        }

        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = engineeringStrain.As3D(mNu);
            break;
        }

        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            engineeringStress.SetZero();
            engineeringStress[0] = (1. - omega) * mE *  elasticEngineeringStrain[0];
            break;
        }

        case NuTo::Constitutive::Output::DAMAGE:
        {
            ConstitutiveIOBase& damage = *itOutput.second;
            damage.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            damage[0] = omega;
            break;
        }

        case NuTo::Constitutive::Output::EXTRAPOLATION_ERROR:
        {
            ConstitutiveIOBase& error = *itOutput.second;
            error.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            error[0] = CalculateStaticDataExtrapolationError(*rElement, rIp, rConstitutiveInput);
            break;
        }

        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "tmp_static_data has to be updated without any other outputs, call it separately.");
            continue;
        }

        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
//        default:
//            throw MechanicsException(__PRETTY_FUNCTION__, "output object " + NuTo::Constitutive::OutputToString(itOutput.first)
//                            + " could not be calculated, check the allocated material law and the section behavior.");
//        }
    }

    //update history variables
    if (performUpdateAtEnd)
        rElement->GetStaticData(rIp)->AsGradientDamage()->SetKappa(currentStaticData.GetKappa());

    return Error::SUCCESSFUL;
}

NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate2D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    // get constitutive inputs
    const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
    const auto& nonlocalEqStrain = *rConstitutiveInput.at(Constitutive::Input::NONLOCAL_EQ_STRAIN);

    auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<2>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

    ConstitutiveStaticDataGradientDamage currentStaticData = GetCurrentStaticData(*rElement, rIp, rConstitutiveInput);
    double omega = CalculateDamage(currentStaticData.GetKappa());

    EquivalentStrainModifiedMises<2> eeq(elasticEngineeringStrain, mCompressiveStrength/mTensileStrength, mNu, rElement->GetSection()->GetType());
    double localEqStrain = eeq.Get();

    bool performUpdateAtEnd = false;


    // calculate coefficients
    double C11, C12, C33;
    switch (rElement->GetSection()->GetType())
    {
    case Section::PLANE_STRAIN:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);
        break;
    case Section::PLANE_STRESS:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients2DPlainStress(mE, mNu);
        break;
    default:
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "[ Invalid type of 2D section behavior found!!!");
    }



    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            engineeringStress[0] = (1 - omega) * (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1]);
            engineeringStress[1] = (1 - omega) * (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0]);
            engineeringStress[2] = (1 - omega) *  C33 * elasticEngineeringStrain[2];
            break;
        }

        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& localEqStrainOut = *itOutput.second;
            localEqStrainOut.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            localEqStrainOut[0] = localEqStrain;
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3,3>(itOutput.first, __PRETTY_FUNCTION__);
            // right coefficients are calculated above
            tangent(0, 0) = (1 - omega) * C11;
            tangent(1, 0) = (1 - omega) * C12;
            tangent(2, 0) = 0;

            tangent(0, 1) = (1 - omega) * C12;
            tangent(1, 1) = (1 - omega) * C11;
            tangent(2, 1) = 0;

            tangent(0, 2) = 0.;
            tangent(1, 2) = 0.;
            tangent(2, 2) = (1 - omega) *C33;
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            if (nonlocalEqStrain[0] == currentStaticData.GetKappa())
            {
                // loading
                double damageDerivative = CalculateDerivativeDamage(nonlocalEqStrain[0]);
                tangent[0] = -damageDerivative * (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1]);
                tangent[1] = -damageDerivative * (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0]);
                tangent[2] = -damageDerivative * (C33 * elasticEngineeringStrain[2]);
            }
            else
            {
                // unloading
                tangent.SetZero();
            }
            break;
        }


        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            tangent = eeq.GetDerivative();
            break;
        }

        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            auto eeqDerivative = eeq.GetDerivative();
            eeqDerivative *= CalculateLocalEqStrainXiFactor(localEqStrain, nonlocalEqStrain[0]);
            tangent = eeqDerivative;
            break;
        }
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveIOBase& xi = *itOutput.second;
            xi.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            xi[0] = CalculateXi(localEqStrain);
            break;
        }

        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
                engineeringStress3D[0] = (1 - omega) * (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1]);
                engineeringStress3D[1] = (1 - omega) * (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0]);
                engineeringStress3D[2] = (1 - omega) *  C12 *(elasticEngineeringStrain[0] + elasticEngineeringStrain[1]);
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = (1 - omega) *  C33 * elasticEngineeringStrain[2];
                break;
            case Section::PLANE_STRESS:
                engineeringStress3D[0] = (1 - omega) * (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1]);
                engineeringStress3D[1] = (1 - omega) * (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0]);
                engineeringStress3D[2] = 0.;
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = (1 - omega) * C33 * elasticEngineeringStrain[2];
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__,"Invalid type of 2D section behavior found!!!");
            }
            break;
        }

        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = engineeringStrain.As3D(mNu, rElement->GetSection()->GetType());
            break;
        }

        case NuTo::Constitutive::Output::DAMAGE:
        {
            ConstitutiveIOBase& damage = *itOutput.second;
            damage.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            damage[0] = omega;
            break;
        }

        case NuTo::Constitutive::Output::EXTRAPOLATION_ERROR:
        {
            ConstitutiveIOBase& error = *itOutput.second;
            error.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            error[0] = CalculateStaticDataExtrapolationError(*rElement, rIp, rConstitutiveInput);
            break;
        }

        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException(__PRETTY_FUNCTION__,"tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            continue;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            continue;
//        default:
//            throw MechanicsException(
//                    std::string(__PRETTY_FUNCTION__,"output object ") + NuTo::Constitutive::OutputToString(itOutput.first)
//                            + std::string(" could not be calculated, check the allocated material law and the section behavior."));
//        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }

    //update history variables
    if (performUpdateAtEnd)
        rElement->GetStaticData(rIp)->AsGradientDamage()->SetKappa(currentStaticData.GetKappa());

    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::GradientDamageEngineeringStress::Evaluate3D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    // get constitutive inputs
    const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain3D();
    const auto& nonlocalEqStrain = *rConstitutiveInput.at(Constitutive::Input::NONLOCAL_EQ_STRAIN);

    auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<3>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

    ConstitutiveStaticDataGradientDamage currentStaticData = GetCurrentStaticData(*rElement, rIp, rConstitutiveInput);
    double omega = CalculateDamage(currentStaticData.GetKappa());

    EquivalentStrainModifiedMises<3> eeq(elasticEngineeringStrain, mCompressiveStrength/mTensileStrength, mNu);
    double localEqStrain = eeq.Get();

    bool performUpdateAtEnd = false;


    double C11, C12, C44;
    std::tie(C11, C12, C44) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            engineeringStress[0] = (1 - omega) * (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[2]);
            engineeringStress[1] = (1 - omega) * (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[2]);
            engineeringStress[2] = (1 - omega) * (C11 * elasticEngineeringStrain[2] + C12 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1]);
            engineeringStress[3] = (1 - omega) *  C44 * elasticEngineeringStrain[3];
            engineeringStress[4] = (1 - omega) *  C44 * elasticEngineeringStrain[4];
            engineeringStress[5] = (1 - omega) *  C44 * elasticEngineeringStrain[5];
            break;
        }

        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = engineeringStrain;
            break;
        }

        case NuTo::Constitutive::Output::LOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& localEqStrainOut = *itOutput.second;
            localEqStrainOut.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            localEqStrainOut[0] = localEqStrain;
            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<6,6>(itOutput.first, __PRETTY_FUNCTION__);

            tangent.SetZero();

            // C11 diagonal:
            tangent(0, 0) = (1 - omega) * C11;
            tangent(1, 1) = (1 - omega) * C11;
            tangent(2, 2) = (1 - omega) * C11;

            // C12 off diagonals:
            tangent(0, 1) = (1 - omega) * C12;
            tangent(0, 2) = (1 - omega) * C12;
            tangent(1, 0) = (1 - omega) * C12;
            tangent(1, 2) = (1 - omega) * C12;
            tangent(2, 0) = (1 - omega) * C12;
            tangent(2, 1) = (1 - omega) * C12;

            // C44 diagonal:
            tangent(3, 3) = (1 - omega) * C44;
            tangent(4, 4) = (1 - omega) * C44;
            tangent(5, 5) = (1 - omega) * C44;
            break;
        }
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            if (nonlocalEqStrain[0] == currentStaticData.GetKappa())
            {
                // loading
                double damageDerivative = CalculateDerivativeDamage(nonlocalEqStrain[0]);
                tangent[0] = -damageDerivative * (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[2] );
                tangent[1] = -damageDerivative * (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[2] );
                tangent[2] = -damageDerivative * (C11 * elasticEngineeringStrain[2] + C12 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1] );

                tangent[3] = -damageDerivative * (C44 * elasticEngineeringStrain[3]);
                tangent[4] = -damageDerivative * (C44 * elasticEngineeringStrain[4]);
                tangent[5] = -damageDerivative * (C44 * elasticEngineeringStrain[5]);

            } else
            {
                // unloading
                tangent.SetZero();
            }
            break;
        }

        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            tangent = eeq.GetDerivative();
            break;
        }

        case NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            auto eeqDerivative = eeq.GetDerivative();
            eeqDerivative *= CalculateLocalEqStrainXiFactor(localEqStrain, nonlocalEqStrain[0]);

            tangent = eeqDerivative;
            break;
        }
        case NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI:
        {
            ConstitutiveIOBase& xi = *itOutput.second;
            xi.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            xi[0] = CalculateXi(localEqStrain);
            break;
        }

        case NuTo::Constitutive::Output::DAMAGE:
        {
            ConstitutiveIOBase& damage = *itOutput.second;
            damage.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            damage[0] = omega;
            break;
        }

        case NuTo::Constitutive::Output::EXTRAPOLATION_ERROR:
        {
            ConstitutiveIOBase& error = *itOutput.second;
            error.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            error[0] = CalculateStaticDataExtrapolationError(*rElement, rIp, rConstitutiveInput);
            break;
        }

        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "tmp_static_data has to be updated without any other outputs, call it separately.");
            continue;
        }

        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
//        default:
//            throw MechanicsException(__PRETTY_FUNCTION__, "output object " + NuTo::Constitutive::OutputToString(itOutput.first)
//                            + " could not be calculated, check the allocated material law and the section behavior.");
//        }
    }

    //update history variables
    if (performUpdateAtEnd)
        rElement->GetStaticData(rIp)->AsGradientDamage()->SetKappa(currentStaticData.GetKappa());

    return Error::SUCCESSFUL;
}

NuTo::ConstitutiveStaticDataGradientDamage NuTo::GradientDamageEngineeringStress::GetCurrentStaticData(ElementBase& rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput) const
{
    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::Input::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "You need to specify the way the static data should be calculated (input list).");

    const auto& calculateStaticData = dynamic_cast<const ConstitutiveCalculateStaticData&>(*itCalculateStaticData->second);

    switch (calculateStaticData.GetCalculateStaticData())
    {
        case CalculateStaticData::USE_PREVIOUS:
        {
            int index = calculateStaticData.GetIndexOfPreviousStaticData();
            return *(rElement.GetStaticDataBase(rIp).GetStaticData(index)->AsGradientDamage());
        }

        case CalculateStaticData::EULER_BACKWARD:
        {
            const auto& nonlocalEqStrain = *rConstitutiveInput.at(Constitutive::Input::NONLOCAL_EQ_STRAIN);

            int index = calculateStaticData.GetIndexOfPreviousStaticData();
            const ConstitutiveStaticDataGradientDamage& oldStaticData = *(rElement.GetStaticDataBase(rIp).GetStaticData(index)->AsGradientDamage());

            ConstitutiveStaticDataGradientDamage newStaticData;
            newStaticData.SetKappa(std::max(nonlocalEqStrain[0], oldStaticData.GetKappa()));

            return newStaticData;
        }

        case CalculateStaticData::EULER_FORWARD:
        {
            auto& staticData = rElement.GetStaticDataBase(rIp);
            assert(staticData.GetNumStaticData() >= 2);

            auto itTimeStep = rConstitutiveInput.find(Constitutive::Input::TIME_STEP);
            if (itTimeStep == rConstitutiveInput.end())
                throw MechanicsException(__PRETTY_FUNCTION__, "TimeStep input needed for EULER_FORWARD.");
            const auto& timeStep = *itTimeStep->second;

            ConstitutiveStaticDataGradientDamage newStaticData;
            double newKappa = ConstitutiveCalculateStaticData::EulerForward(
                    staticData.GetStaticData(1)->AsGradientDamage()->GetKappa(),
                    staticData.GetStaticData(2)->AsGradientDamage()->GetKappa(),
                    timeStep);

//            std::cout << newKappa << std::endl;

            newStaticData.SetKappa(newKappa);
            return newStaticData;
        }

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Cannot calculate the static data in the requested way.");
    }

}

double NuTo::GradientDamageEngineeringStress::CalculateStaticDataExtrapolationError(ElementBase& rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput) const
{
    // static data 0 contains the extrapolated values \tilde \kappa_n
    // static data 1 contains the implicit data \kappa_n-1
    // static data 2 contains the implicit data \kappa_n-2

    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::Input::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "You need to specify the way the static data should be calculated (input list).");

    const auto& ipData = rElement.GetStaticDataBase(rIp);

    double eeq = (*rConstitutiveInput.at(Constitutive::Input::NONLOCAL_EQ_STRAIN))[0];
    double k_n_t = ipData.GetStaticData(0)->AsGradientDamage()->GetKappa();
    double k_n_1 = ipData.GetStaticData(1)->AsGradientDamage()->GetKappa();

    double k_n = std::max(eeq, k_n_1); // calculate kappa implicitly
    return std::abs(k_n - k_n_t);
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticData1D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage;
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticData2D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage;
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::GradientDamageEngineeringStress::AllocateStaticData3D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataGradientDamage;
}


//! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
//! @param rDofRow ... row dof
//! @param rDofCol ... column dof
//! @param rTimeDerivative ... time derivative
bool NuTo::GradientDamageEngineeringStress::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    assert(rTimeDerivative>-1);
    if (rTimeDerivative<1)
    {
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::DISPLACEMENTS,      Node::DISPLACEMENTS):
        case Node::CombineDofs(Node::DISPLACEMENTS,      Node::NONLOCALEQSTRAIN):
        case Node::CombineDofs(Node::NONLOCALEQSTRAIN,   Node::DISPLACEMENTS):
        case Node::CombineDofs(Node::NONLOCALEQSTRAIN,   Node::NONLOCALEQSTRAIN):
            return true;
        default:
            return false;
        }
    }
    return false;
}

// parameters /////////////////////////////////////////////////////////////

double NuTo::GradientDamageEngineeringStress::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        return mCompressiveStrength;
    case Constitutive::eConstitutiveParameter::DENSITY:
        return mRho;
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        return mFractureEnergy;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        return mNonlocalRadius;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
        return mNonlocalRadiusParameter;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        return mNu;
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
        return mTensileStrength;
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        return mThermalExpansionCoefficient;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        return mE;
    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
        return mDamageLawType;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Constitutive law does not have the requested variable");
    }
}

void NuTo::GradientDamageEngineeringStress::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        mCompressiveStrength = rValue;
        break;
    case Constitutive::eConstitutiveParameter::DENSITY:
        mRho = rValue;
        break;
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        mFractureEnergy = rValue;
        break;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        mNonlocalRadius = rValue;
        break;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER:
        mNonlocalRadiusParameter = rValue;
        break;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        mNu = rValue;
        break;
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
        mTensileStrength = rValue;
        break;
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        mThermalExpansionCoefficient = rValue;
        break;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        mE = rValue;
        break;
    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
        mDamageLawType = (Constitutive::eDamageLawType) rValue;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
    SetParametersValid();
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
    case NuTo::Element::CONTINUUMELEMENT:
        return true;
    case NuTo::Element::CONTINUUMBOUNDARYELEMENT:
        return true;

    default:
        return false;
    }
}




//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
//! @param rLogger stream for the output
void NuTo::GradientDamageEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus         : " << mE << "\n";
    rLogger << "    Poisson's ratio         : " << mNu << "\n";
    rLogger << "    thermal expansion coeff : " << mThermalExpansionCoefficient << "\n";
    rLogger << "    nonlocal radius         : " << mNonlocalRadius << "\n";
    rLogger << "    nonlocal radius param   : " << mNonlocalRadiusParameter << "\n";
    rLogger << "    tensile strength        : " << mTensileStrength << "\n";
    rLogger << "    compressive strength    : " << mCompressiveStrength << "\n";
    rLogger << "    fracture energy         : " << mFractureEnergy << "\n";
}

// check parameters
void NuTo::GradientDamageEngineeringStress::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::DENSITY, mRho);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, mE);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, mNu);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, mNonlocalRadius);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, mNonlocalRadiusParameter);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, mTensileStrength);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, mFractureEnergy);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, mThermalExpansionCoefficient);
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


double NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainXiFactor(double rLocalEqStrain, double rNonlocalEqStrain) const
{
    double xi = mNonlocalRadius;
    double dXi = 0;

    if (mNonlocalRadiusParameter != 0.)
    {
        double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
        if (rLocalEqStrain <= e_xi)
        {
            double c0 = mNonlocalRadius/100.;
            xi = c0 + (mNonlocalRadius - c0) * (rLocalEqStrain / e_xi);
            dXi = (mNonlocalRadius - c0) / e_xi;
        }
    }
    return 1. / xi - (rLocalEqStrain - rNonlocalEqStrain) / (xi * xi) * dXi;
}

double NuTo::GradientDamageEngineeringStress::CalculateXi(double rLocalEqStrain) const
{
    double xi = mNonlocalRadius;
    if (mNonlocalRadiusParameter != 0.)
    {
        double e_xi = mTensileStrength / mE * mNonlocalRadiusParameter;
        if (rLocalEqStrain <= e_xi)
        {
            double c0 = mNonlocalRadius/100.;
            xi = c0 + (mNonlocalRadius - c0) * (rLocalEqStrain / e_xi);
        }
    }
    return xi;
}
