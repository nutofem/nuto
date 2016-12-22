#include "mechanics/constitutive/laws/LocalDamageModel.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "base/ErrorEnum.h"
#include "base/Logger.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/sections/SectionBase.h"
#include "mechanics/sections/SectionEnum.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveMatrixXd.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "mechanics/constitutive/inputoutput/EquivalentStrain.h"


using NuTo::Constitutive::eDamageLawType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Constitutive::eOutput;
using NuTo::Constitutive::eInput;
using NuTo::Constitutive::eConstitutiveType;


NuTo::LocalDamageModel::LocalDamageModel() :
    ConstitutiveBase(),
    mDensity(0.),
    mYoungsModulus(0.),
    mPoissonsRatio(0.),
    mTensileStrength(0.),
    mCompressiveStrength(0.),
    mFractureEnergy(0.),
    mDamageLawType(eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD)
{
}


double NuTo::LocalDamageModel::GetParameterDouble(eConstitutiveParameter rIdentifier) const
{
    CheckParameters();

    switch (rIdentifier)
    {
    case eConstitutiveParameter::DENSITY:
        return mDensity;

    case eConstitutiveParameter::YOUNGS_MODULUS:
        return mYoungsModulus;

    case eConstitutiveParameter::POISSONS_RATIO:
        return mPoissonsRatio;

    case eConstitutiveParameter::TENSILE_STRENGTH:
        return mTensileStrength;

    case eConstitutiveParameter::FRACTURE_ENERGY:
        return mFractureEnergy;

    case eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        return mCompressiveStrength;

    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
        return static_cast<double>(mDamageLawType);

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Constitutive law does not have the requested variable");
    }
}

void NuTo::LocalDamageModel::SetParameterDouble(eConstitutiveParameter rIdentifier, double rValue)
{
    switch (rIdentifier)
    {
    case eConstitutiveParameter::COMPRESSIVE_STRENGTH:
        mCompressiveStrength = rValue;
        break;

    case eConstitutiveParameter::FRACTURE_ENERGY:
        mFractureEnergy = rValue;
        break;

    case eConstitutiveParameter::DENSITY:
        mDensity = rValue;
        break;

    case eConstitutiveParameter::YOUNGS_MODULUS:
        mYoungsModulus = rValue;
        break;

    case eConstitutiveParameter::POISSONS_RATIO:
        mPoissonsRatio = rValue;
        break;

    case eConstitutiveParameter::TENSILE_STRENGTH:
        mTensileStrength = rValue;
        break;

    case Constitutive::eConstitutiveParameter::DAMAGE_LAW:
        mDamageLawType = (Constitutive::eDamageLawType) rValue;
        break;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Constitutive law does not have the requested variable");
    }

}

//! @brief ... get type of constitutive relationship
eConstitutiveType NuTo::LocalDamageModel::GetType() const
{
    return eConstitutiveType::LOCAL_DAMAGE_MODEL;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
bool NuTo::LocalDamageModel::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
        return true;
    default:
        return false;
    }
}



//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::LocalDamageModel::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "Info function not yet implemented" << "\n";

}


void NuTo::LocalDamageModel::CheckParameters() const {}


double NuTo::LocalDamageModel::Evaluate2D(NuTo::ConstitutivePlaneState planeState, double oldKappa,
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    // get constitutive inputs
    const auto& elasticEngineeringStrain = rConstitutiveInput.at(eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();

    EquivalentStrainModifiedMises<2> eeq(elasticEngineeringStrain, mCompressiveStrength/mTensileStrength,
            mPoissonsRatio, planeState.GetPlaneState());
    double localEqStrain = eeq.Get();

    double kappa = std::max(localEqStrain, oldKappa);
    double omega = CalculateDamage(kappa);

    bool performUpdateAtEnd = false;

    // calculate coefficients
    double C11, C12, C33;
    switch (planeState.GetPlaneState())
    {
    case ePlaneState::PLANE_STRAIN:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mYoungsModulus, mPoissonsRatio);
        break;
    case ePlaneState::PLANE_STRESS:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients2DPlaneStress(mYoungsModulus, mPoissonsRatio);
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Invalid type of 2D section behavior found.");
    }



    /////////////////////////////////////////////////
    //         LOOP OVER OUTPUT REQUESTS           //
    /////////////////////////////////////////////////

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {

        case eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            engineeringStress[0] = (1 - omega) * (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1]);
            engineeringStress[1] = (1 - omega) * (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0]);
            engineeringStress[2] = (1 - omega) *  C33 * elasticEngineeringStrain[2];
            break;
        }

        case eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3,3>(itOutput.first, __PRETTY_FUNCTION__);

            double  dDamageDLocalEqStrain = CalculateDerivativeDamage(kappa);
            auto    dLocalEqStrainDStrain = eeq.GetDerivative();

            std::vector<double> effectiveStress(3);
            effectiveStress[0] = (C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1]);
            effectiveStress[1] = (C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0]);
            effectiveStress[2] =  C33 * elasticEngineeringStrain[2];

            // right coefficients are calculated above
            tangent(0, 0) = (1 - omega) * C11 - dDamageDLocalEqStrain * dLocalEqStrainDStrain[0] * effectiveStress[0];
            tangent(1, 0) = (1 - omega) * C12 - dDamageDLocalEqStrain * dLocalEqStrainDStrain[1] * effectiveStress[0];
            tangent(2, 0) = 0.                - dDamageDLocalEqStrain * dLocalEqStrainDStrain[2] * effectiveStress[0];

            tangent(0, 1) = (1 - omega) * C12 - dDamageDLocalEqStrain * dLocalEqStrainDStrain[0] * effectiveStress[1];
            tangent(1, 1) = (1 - omega) * C11 - dDamageDLocalEqStrain * dLocalEqStrainDStrain[1] * effectiveStress[1];
            tangent(2, 1) = 0.                - dDamageDLocalEqStrain * dLocalEqStrainDStrain[2] * effectiveStress[1];

            tangent(0, 2) = 0.                - dDamageDLocalEqStrain * dLocalEqStrainDStrain[0] * effectiveStress[2];
            tangent(1, 2) = 0.                - dDamageDLocalEqStrain * dLocalEqStrainDStrain[1] * effectiveStress[2];
            tangent(2, 2) = (1 - omega) *C33  - dDamageDLocalEqStrain * dLocalEqStrainDStrain[2] * effectiveStress[2];
            break;
        }

        case eOutput::DAMAGE:
        {
            ConstitutiveIOBase& damage = *itOutput.second;
            damage.AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            damage[0] = omega;
            break;
        }

        case eOutput::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException(__PRETTY_FUNCTION__,"tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            continue;
        case eOutput::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            continue;

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }


    // return old/new history variables
    if (performUpdateAtEnd)
    {
        return kappa;
    }
    else
    {
        return oldKappa;
    }
}

template <int TDim>
NuTo::eError NuTo::LocalDamageModel::Evaluate(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput,
    Data& rStaticData)
{
    if (TDim != 2)
        throw MechanicsException(__PRETTY_FUNCTION__, "IMPLEMENT ME!!!");

    auto itCalculateStaticData = rConstitutiveInput.find(eInput::CALCULATE_STATIC_DATA);
    const auto
        & calculateStaticData = dynamic_cast<const ConstitutiveCalculateStaticData&>(*itCalculateStaticData->second);
    int index = calculateStaticData.GetIndexOfPreviousStaticData();

    double oldKappa = rStaticData.GetData(index);

    const auto& planeState =
        *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());
    double newKappa = Evaluate2D(planeState, oldKappa, rConstitutiveInput, rConstitutiveOutput);

    // update history variables
    rStaticData.GetData() = newKappa;
    return eError::SUCCESSFUL;
}

NuTo::ConstitutiveInputMap NuTo::LocalDamageModel::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    constitutiveInputMap[eInput::ENGINEERING_STRAIN];

    return constitutiveInputMap;
}


double NuTo::LocalDamageModel::GetCurrentStaticData(Data& rStaticData,
        const ConstitutiveInputMap& rConstitutiveInput) const
{
    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::eInput::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw MechanicsException(__PRETTY_FUNCTION__,
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
                throw MechanicsException(__PRETTY_FUNCTION__, "TimeStep input needed for EULER_FORWARD.");
            const auto& timeStep = *(itTimeStep->second);

            assert(rStaticData.GetNumData() >= 2);

            return ConstitutiveCalculateStaticData::EulerForward(
                    rStaticData.GetData(1), rStaticData.GetData(2), timeStep);
        }

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Cannot calculate the static data in the requested way.");
    }

}


double NuTo::LocalDamageModel::CalculateDamage(double rKappa) const
{
    double omega = 0;
    double e_0 = mTensileStrength / mYoungsModulus;
    double e_f = mFractureEnergy / mTensileStrength;

    switch (mDamageLawType)
    {
    case eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
    {
        if (rKappa > e_0)
        {
            omega = 1 - (1-mMaxDamage) - mMaxDamage*e_0 / rKappa * (1 - mAlpha + mAlpha * std::exp((e_0 - rKappa) / e_f));
        }
        break;
    }
    case Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING:
    {
        if (rKappa > e_0)
        {
            omega = 1 - e_0 / rKappa * std::exp((e_0 - rKappa) / e_f);
        }
        break;
    }

    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "The required damage law is not implemented. ");
    }

    return omega;
}

double NuTo::LocalDamageModel::CalculateDerivativeDamage(double rKappa) const
{
    double DomegaDkappa = 0;

    double e_0 = mTensileStrength / mYoungsModulus;
    double e_f = mFractureEnergy / mTensileStrength;

    switch (mDamageLawType)
    {
    case eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD:
    {
        if (rKappa > e_0)
        {
            DomegaDkappa = mMaxDamage * e_0 / rKappa
                * ((1 / rKappa + 1 / e_f) * mAlpha * std::exp((e_0 - rKappa) / e_f) + (1 - mAlpha) / rKappa);
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

    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "The required damage law is not implemented. ");
    }

    return DomegaDkappa;
}


#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::LocalDamageModel::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::LocalDamageModel::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::LocalDamageModel::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::LocalDamageModel::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::LocalDamageModel::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::LocalDamageModel::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::LocalDamageModel::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LocalDamageModel" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
            & BOOST_SERIALIZATION_NVP(mMaxBondStress)
            & BOOST_SERIALIZATION_NVP(mResidualBondStress)
            & BOOST_SERIALIZATION_NVP(mSlipAtMaxBondStress)
            & BOOST_SERIALIZATION_NVP(mSlipAtResidualBondStress)
            & BOOST_SERIALIZATION_NVP(mNormalStiffness)
            & BOOST_SERIALIZATION_NVP(mAlpha);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LocalDamageModel" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LocalDamageModel)
#endif // ENABLE_SERIALIZATION


bool NuTo::LocalDamageModel::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    return rTimeDerivative<1 &&
        rDofRow == Node::eDof::DISPLACEMENTS &&
        rDofCol ==Node::eDof::DISPLACEMENTS;
}

template NuTo::eError NuTo::LocalDamageModel::Evaluate<1>(
    const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput, Data& rStaticData);
template NuTo::eError NuTo::LocalDamageModel::Evaluate<2>(
    const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput, Data& rStaticData);
template NuTo::eError NuTo::LocalDamageModel::Evaluate<3>(
    const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput, Data& rStaticData);