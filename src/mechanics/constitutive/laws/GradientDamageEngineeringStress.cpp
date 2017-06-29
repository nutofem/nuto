#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/EquivalentStrain.h"

#include "base/Logger.h"
#include "base/Exception.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/constitutive/damageLaws/DamageLaw.h"

#include "mechanics/timeIntegration/ImplExCallback.h"

NuTo::GradientDamageEngineeringStress::GradientDamageEngineeringStress()
    : ConstitutiveBase()
    , mRho(0.)
    , mE(0.)
    , mNu(0.)
    , mNonlocalRadius(0.)
    , mThermalExpansionCoefficient(0.)
    , mTensileStrength(0.)
    , mCompressiveStrength(0.)
    , mDamageLaw()
    , mImplExCallback(std::make_shared<ImplExCallback>())
{
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::binary_oarchive& ar,
                                                               const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::binary_iarchive& ar,
                                                               const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::xml_oarchive& ar,
                                                               const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::xml_iarchive& ar,
                                                               const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::text_oarchive& ar,
                                                               const unsigned int version);
template void NuTo::GradientDamageEngineeringStress::serialize(boost::archive::text_iarchive& ar,
                                                               const unsigned int version);
template <class Archive>
void NuTo::GradientDamageEngineeringStress::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize GradientDamageEngineeringStress"
              << "\n";
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase) & BOOST_SERIALIZATION_NVP(mRho) &
            BOOST_SERIALIZATION_NVP(mE) & BOOST_SERIALIZATION_NVP(mNu) & BOOST_SERIALIZATION_NVP(mNonlocalRadius) &
            BOOST_SERIALIZATION_NVP(mNonlocalRadiusParameter) & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient) &
            BOOST_SERIALIZATION_NVP(mTensileStrength) & BOOST_SERIALIZATION_NVP(mCompressiveStrength) &
            BOOST_SERIALIZATION_NVP(mFractureEnergy) & BOOST_SERIALIZATION_NVP(mDamageLawType);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GradientDamageEngineeringStress \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::GradientDamageEngineeringStress)
#endif // ENABLE_SERIALIZATION

NuTo::ConstitutiveInputMap
NuTo::GradientDamageEngineeringStress::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                             const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    // Always allocate these two guys. For a dof separated approach, this could be handled more efficiently.
    constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
    constitutiveInputMap[Constitutive::eInput::NONLOCAL_EQ_STRAIN];

    return constitutiveInputMap;
}

template double
NuTo::GradientDamageEngineeringStress::EvaluateStaticData<1>(const ConstitutiveInputMap& rConstitutiveInput,
                                                             const ConstitutiveOutputMap& rConstitutiveOutput,
                                                             Data& rStaticData);
template double
NuTo::GradientDamageEngineeringStress::EvaluateStaticData<2>(const ConstitutiveInputMap& rConstitutiveInput,
                                                             const ConstitutiveOutputMap& rConstitutiveOutput,
                                                             Data& rStaticData);
template double
NuTo::GradientDamageEngineeringStress::EvaluateStaticData<3>(const ConstitutiveInputMap& rConstitutiveInput,
                                                             const ConstitutiveOutputMap& rConstitutiveOutput,
                                                             Data& rStaticData);

template <int TDim>
double NuTo::GradientDamageEngineeringStress::EvaluateStaticData(const ConstitutiveInputMap& rConstitutiveInput,
                                                                 const ConstitutiveOutputMap& rConstitutiveOutput,
                                                                 Data& rStaticData)
{
    double kappa = GetCurrentStaticData(rStaticData, rConstitutiveInput);
    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::EXTRAPOLATION_ERROR:
        {
            ConstitutiveIOBase& error = *itOutput.second;
            error.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            error[0] = CalculateStaticDataExtrapolationError(rStaticData, rConstitutiveInput);
            error.SetIsCalculated(true);
            break;
        }

        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        {
            throw Exception(
                    __PRETTY_FUNCTION__,
                    "tmp_static_data has to be updated without any other outputs, call it separately.");
        }

        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            rStaticData.SetData(kappa);
        }
        default:
            continue;
        }
    }

    return kappa;
}


namespace NuTo // template specialization in same namespace as definition
{
template <>
void NuTo::GradientDamageEngineeringStress::EvaluateWithKappa<1>(const ConstitutiveInputMap& rConstitutiveInput,
                                                                 const ConstitutiveOutputMap& rConstitutiveOutput,
                                                                 StaticDataType rKappa, double rKappaTangent)
{
    // get constitutive inputs
    const auto& engineeringStrain =
            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D();

    double omega = mDamageLaw->CalculateDamage(rKappa);

    EquivalentStrainModifiedMises<1> eeq(engineeringStrain, mCompressiveStrength / mTensileStrength, mNu);
    double localEqStrain = eeq.Get();

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            engineeringStress[0] = (1. - omega) * mE * engineeringStrain[0];
            break;
        }

        case NuTo::Constitutive::eOutput::LOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& localEqStrainOut = *itOutput.second;
            localEqStrainOut.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            localEqStrainOut[0] = localEqStrain;
            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1, 1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = (1. - omega) * mE;
            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);
            if (rKappaTangent == 0.)
            {
                tangent.SetZero();
                break;
            }
            tangent(0, 0) = -rKappaTangent * mDamageLaw->CalculateDerivative(rKappa) * mE * engineeringStrain[0];
            break;
        }

        case NuTo::Constitutive::eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);
            tangent[0] = eeq.GetDerivative()[0];
            break;
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = engineeringStrain.As3D(mNu);
            break;
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            engineeringStress.SetZero();
            engineeringStress[0] = (1. - omega) * mE * engineeringStrain[0];
            break;
        }

        case NuTo::Constitutive::eOutput::DAMAGE:
        {
            ConstitutiveIOBase& damage = *itOutput.second;
            damage.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            damage[0] = omega;
            break;
        }

        case NuTo::Constitutive::eOutput::NONLOCAL_RADIUS:
        {
            ConstitutiveIOBase& radius = *itOutput.second;
            radius.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            radius[0] = mNonlocalRadius;
            break;
        }

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
        //        default:
        //            throw Exception(__PRETTY_FUNCTION__, "output object " +
        //            NuTo::Constitutive::OutputToString(itOutput.first)
        //                            + " could not be calculated, check the allocated material law and the section
        //                            behavior.");
        //        }
    }
}

template <>
void NuTo::GradientDamageEngineeringStress::EvaluateWithKappa<2>(const ConstitutiveInputMap& rConstitutiveInput,
                                                                 const ConstitutiveOutputMap& rConstitutiveOutput,
                                                                 StaticDataType rKappa, double rKappaTangent)
{
    // get constitutive inputs
    const auto& engineeringStrain =
            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
    const auto& planeState =
            *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());

    double omega = mDamageLaw->CalculateDamage(rKappa);

    EquivalentStrainModifiedMises<2> eeq(engineeringStrain, mCompressiveStrength / mTensileStrength, mNu,
                                         planeState.GetPlaneState());
    double localEqStrain = eeq.Get();

    // calculate coefficients
    double C11, C12, C33;
    switch (planeState.GetPlaneState())
    {
    case ePlaneState::PLANE_STRAIN:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);
        break;
    case ePlaneState::PLANE_STRESS:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients2DPlaneStress(mE, mNu);
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Invalid type of 2D section behavior found.");
    }


    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            engineeringStress[0] = (1 - omega) * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
            engineeringStress[1] = (1 - omega) * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
            engineeringStress[2] = (1 - omega) * C33 * engineeringStrain[2];
            break;
        }

        case NuTo::Constitutive::eOutput::LOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& localEqStrainOut = *itOutput.second;
            localEqStrainOut.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            localEqStrainOut[0] = localEqStrain;
            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);
            // right coefficients are calculated above
            tangent(0, 0) = (1 - omega) * C11;
            tangent(1, 0) = (1 - omega) * C12;
            tangent(2, 0) = 0;

            tangent(0, 1) = (1 - omega) * C12;
            tangent(1, 1) = (1 - omega) * C11;
            tangent(2, 1) = 0;

            tangent(0, 2) = 0.;
            tangent(1, 2) = 0.;
            tangent(2, 2) = (1 - omega) * C33;
            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            if (rKappaTangent == 0.)
            {
                tangent.SetZero();
                break;
            }
            double damageDerivative = rKappaTangent * mDamageLaw->CalculateDerivative(rKappa);
            tangent[0] = -damageDerivative * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
            tangent[1] = -damageDerivative * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
            tangent[2] = -damageDerivative * (C33 * engineeringStrain[2]);
            break;
        }


        case NuTo::Constitutive::eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            tangent = eeq.GetDerivative();
            break;
        }


        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            switch (planeState.GetPlaneState())
            {
            case ePlaneState::PLANE_STRAIN:
                engineeringStress3D[0] = (1 - omega) * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
                engineeringStress3D[1] = (1 - omega) * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
                engineeringStress3D[2] = (1 - omega) * C12 * (engineeringStrain[0] + engineeringStrain[1]);
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = (1 - omega) * C33 * engineeringStrain[2];
                break;
            case ePlaneState::PLANE_STRESS:
                engineeringStress3D[0] = (1 - omega) * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
                engineeringStress3D[1] = (1 - omega) * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
                engineeringStress3D[2] = 0.;
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = (1 - omega) * C33 * engineeringStrain[2];
                break;
            default:
                throw Exception(__PRETTY_FUNCTION__, "Invalid type of 2D section behavior found!!!");
            }
            break;
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = engineeringStrain.As3D(mNu, planeState.GetPlaneState());
            break;
        }

        case NuTo::Constitutive::eOutput::DAMAGE:
        {
            ConstitutiveIOBase& damage = *itOutput.second;
            damage.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            damage[0] = omega;
            break;
        }

        case NuTo::Constitutive::eOutput::NONLOCAL_RADIUS:
        {
            ConstitutiveIOBase& radius = *itOutput.second;
            radius.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            radius[0] = mNonlocalRadius;
            break;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

template <>
void NuTo::GradientDamageEngineeringStress::EvaluateWithKappa<3>(const ConstitutiveInputMap& rConstitutiveInput,
                                                                 const ConstitutiveOutputMap& rConstitutiveOutput,
                                                                 StaticDataType rKappa, double rKappaTangent)
{
    // get constitutive inputs
    const auto& engineeringStrain =
            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain3D();

    double omega = mDamageLaw->CalculateDamage(rKappa);

    EquivalentStrainModifiedMises<3> eeq(engineeringStrain, mCompressiveStrength / mTensileStrength, mNu);
    double localEqStrain = eeq.Get();

    double C11, C12, C44;
    std::tie(C11, C12, C44) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            engineeringStress[0] = (1 - omega) * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1] +
                                                  C12 * engineeringStrain[2]);
            engineeringStress[1] = (1 - omega) * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0] +
                                                  C12 * engineeringStrain[2]);
            engineeringStress[2] = (1 - omega) * (C11 * engineeringStrain[2] + C12 * engineeringStrain[0] +
                                                  C12 * engineeringStrain[1]);
            engineeringStress[3] = (1 - omega) * C44 * engineeringStrain[3];
            engineeringStress[4] = (1 - omega) * C44 * engineeringStrain[4];
            engineeringStress[5] = (1 - omega) * C44 * engineeringStrain[5];
            break;
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = engineeringStrain;
            break;
        }

        case NuTo::Constitutive::eOutput::LOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& localEqStrainOut = *itOutput.second;
            localEqStrainOut.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            localEqStrainOut[0] = localEqStrain;
            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<6, 6>(itOutput.first, __PRETTY_FUNCTION__);

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
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            if (rKappaTangent == 0.)
            {
                tangent.SetZero();
                break;
            }
            double damageDerivative = rKappaTangent * mDamageLaw->CalculateDerivative(rKappa);
            tangent[0] = -damageDerivative *
                         (C11 * engineeringStrain[0] + C12 * engineeringStrain[1] + C12 * engineeringStrain[2]);
            tangent[1] = -damageDerivative *
                         (C11 * engineeringStrain[1] + C12 * engineeringStrain[0] + C12 * engineeringStrain[2]);
            tangent[2] = -damageDerivative *
                         (C11 * engineeringStrain[2] + C12 * engineeringStrain[0] + C12 * engineeringStrain[1]);

            tangent[3] = -damageDerivative * (C44 * engineeringStrain[3]);
            tangent[4] = -damageDerivative * (C44 * engineeringStrain[4]);
            tangent[5] = -damageDerivative * (C44 * engineeringStrain[5]);
            break;
        }

        case NuTo::Constitutive::eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            tangent = eeq.GetDerivative();
            break;
        }


        case NuTo::Constitutive::eOutput::DAMAGE:
        {
            ConstitutiveIOBase& damage = *itOutput.second;
            damage.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            damage[0] = omega;
            break;
        }

        case NuTo::Constitutive::eOutput::NONLOCAL_RADIUS:
        {
            ConstitutiveIOBase& radius = *itOutput.second;
            radius.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            radius[0] = mNonlocalRadius;
            break;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
        //        default:
        //            throw Exception(__PRETTY_FUNCTION__, "output object " +
        //            NuTo::Constitutive::OutputToString(itOutput.first)
        //                            + " could not be calculated, check the allocated material law and the section
        //                            behavior.");
        //        }
    }
}

} // namespace NuTo

double NuTo::GradientDamageEngineeringStress::GetCurrentStaticData(Data& rStaticData,
                                                                   const ConstitutiveInputMap& rConstitutiveInput) const
{
    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::eInput::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw Exception(__PRETTY_FUNCTION__,
                                 "You need to specify the way the static data should be calculated (input list).");

    const auto& calculateStaticData =
            *static_cast<const ConstitutiveCalculateStaticData*>(itCalculateStaticData->second.get());

    switch (calculateStaticData.GetCalculateStaticData())
    {
    case eCalculateStaticData::USE_PREVIOUS:
    {
        int index = calculateStaticData.GetIndexOfPreviousStaticData();
        return rStaticData.GetData(index);
    }
    case eCalculateStaticData::EULER_BACKWARD:
    {

        int index = calculateStaticData.GetIndexOfPreviousStaticData();
        double oldKappa = rStaticData.GetData(index);
        const auto& nonlocalEqStrain = *rConstitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN);
        return std::max(nonlocalEqStrain[0], oldKappa);
    }

    case eCalculateStaticData::EULER_FORWARD:
    {
        auto itTimeStep = rConstitutiveInput.find(Constitutive::eInput::TIME_STEP);
        if (itTimeStep == rConstitutiveInput.end())
            throw Exception(__PRETTY_FUNCTION__, "TimeStep input needed for EULER_FORWARD.");
        const auto& timeStep = *(itTimeStep->second);

        assert(rStaticData.GetNumData() >= 2);

        return ConstitutiveCalculateStaticData::EulerForward(rStaticData.GetData(1), rStaticData.GetData(2), timeStep);
    }

    default:
        throw Exception(__PRETTY_FUNCTION__, "Cannot calculate the static data in the requested way.");
    }
}

double NuTo::GradientDamageEngineeringStress::CalculateStaticDataExtrapolationError(
        Data& rStaticData, const ConstitutiveInputMap& rConstitutiveInput) const
{
    // static data 0 contains the extrapolated values \tilde \kappa_n
    // static data 1 contains the implicit data \kappa_n-1
    // static data 2 contains the implicit data \kappa_n-2

    double eeq = (*rConstitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN))[0];
    double k_n_t = rStaticData.GetData(0);
    double k_nn = rStaticData.GetData(1);
    double k_n = std::max(eeq, k_nn); // calculate kappa implicitly

    return mImplExCallback->GetError(k_n_t, k_n, k_nn);
}


bool NuTo::GradientDamageEngineeringStress::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow,
                                                                          NuTo::Node::eDof rDofCol,
                                                                          int rTimeDerivative) const
{
    assert(rTimeDerivative > -1);
    if (rTimeDerivative < 1)
    {
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN):
        case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::NONLOCALEQSTRAIN):
            return true;
        default:
            return false;
        }
    }
    return false;
}

// parameters /////////////////////////////////////////////////////////////

double
NuTo::GradientDamageEngineeringStress::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        return mCompressiveStrength;
    case Constitutive::eConstitutiveParameter::DENSITY:
        return mRho;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        return mNonlocalRadius;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        return mNu;
    case Constitutive::eConstitutiveParameter::TENSILE_STRENGTH:
        return mTensileStrength;
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        return mThermalExpansionCoefficient;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        return mE;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void NuTo::GradientDamageEngineeringStress::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier,
                                                               double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        mCompressiveStrength = rValue;
        break;
    case Constitutive::eConstitutiveParameter::DENSITY:
        mRho = rValue;
        break;
    case Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS:
        mNonlocalRadius = rValue;
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
    default:
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
    SetParametersValid();
}

///////////////////////////////////////////////////////////////////////////

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::GradientDamageEngineeringStress::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS;
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
    rLogger << "    tensile strength        : " << mTensileStrength << "\n";
    rLogger << "    compressive strength    : " << mCompressiveStrength << "\n";
}

// check parameters
void NuTo::GradientDamageEngineeringStress::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::DENSITY, mRho);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, mE);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, mNu);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, mNonlocalRadius);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, mTensileStrength);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT,
                                           mThermalExpansionCoefficient);
}
template void
NuTo::GradientDamageEngineeringStress::Evaluate<1>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                   const NuTo::ConstitutiveOutputMap& rConstitutiveOutput,
                                                   NuTo::GradientDamageEngineeringStress::Data& rStaticData);
template void
NuTo::GradientDamageEngineeringStress::Evaluate<2>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                   const NuTo::ConstitutiveOutputMap& rConstitutiveOutput,
                                                   NuTo::GradientDamageEngineeringStress::Data& rStaticData);
template void
NuTo::GradientDamageEngineeringStress::Evaluate<3>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                   const NuTo::ConstitutiveOutputMap& rConstitutiveOutput,
                                                   NuTo::GradientDamageEngineeringStress::Data& rStaticData);

template <int TDim>
void NuTo::GradientDamageEngineeringStress::Evaluate(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                     const NuTo::ConstitutiveOutputMap& rConstitutiveOutput,
                                                     NuTo::GradientDamageEngineeringStress::Data& rStaticData)
{
    // this split allows reusing the EvaluteWithKappa from other classes
    double kappa = EvaluateStaticData<TDim>(rConstitutiveInput, rConstitutiveOutput, rStaticData);
    double nonlocalEqStrain = rConstitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN)->operator[](0);
    double kappaTangent = kappa == nonlocalEqStrain ? 1.0 : 0.0; // = 1 true, 0 false. perfect tangent.
    EvaluateWithKappa<TDim>(rConstitutiveInput, rConstitutiveOutput, kappa, kappaTangent);
}

void NuTo::GradientDamageEngineeringStress::SetExtrapolation(std::shared_ptr<ImplExCallback> rCallback)
{
    mImplExCallback = rCallback;
}
