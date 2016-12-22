#include "mechanics/constitutive/laws/FibreMatrixBondStressSlip.h"

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

#include "mechanics/constitutive/inputoutput/ConstitutiveMatrixXd.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"

NuTo::FibreMatrixBondStressSlip::FibreMatrixBondStressSlip(int dimension) :
        ConstitutiveBase(),
        mMaxBondStress(0),
        mResidualBondStress(0),
        mSlipAtMaxBondStress(0),
        mSlipAtResidualBondStress(0),
        mAlpha(0),
        mNormalStiffness(0),
        mGlobalDimension(dimension)
{
}




//! @brief ... gets a variable of the constitutive law which is selected by an enum
double NuTo::FibreMatrixBondStressSlip::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    CheckParameters();

    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS:
        return this->mNormalStiffness;
    case Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS:
        return this->mResidualBondStress;
    case Constitutive::eConstitutiveParameter::MAX_BOND_STRESS:
        return this->mMaxBondStress;
    case Constitutive::eConstitutiveParameter::ALPHA:
        return this->mAlpha;
    case Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS:
        return this->mSlipAtMaxBondStress;
    case Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS:
        return this->mSlipAtResidualBondStress;
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive law does not have the requested variable");
    }
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
void NuTo::FibreMatrixBondStressSlip::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS:
    {
        mNormalStiffness = rValue;
        break;
    }
    case Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS:
    {
        mResidualBondStress = rValue;
        break;
    }
    case Constitutive::eConstitutiveParameter::MAX_BOND_STRESS:
    {
        mMaxBondStress = rValue;
        break;
    }
    case Constitutive::eConstitutiveParameter::ALPHA:
    {
        mAlpha = rValue;
        break;
    }
    case Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS:
    {
        mSlipAtMaxBondStress = rValue;
        break;
    }
    case Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS:
    {
        mSlipAtResidualBondStress = rValue;
        break;
    }
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive law does not have the requested variable");
    }

}

//! @brief ... get type of constitutive relationship
NuTo::Constitutive::eConstitutiveType NuTo::FibreMatrixBondStressSlip::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
bool NuTo::FibreMatrixBondStressSlip::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::ELEMENT2DINTERFACE:
        return true;
    default:
        return false;
    }
}



//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::FibreMatrixBondStressSlip::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Normal stiffness               : " << mNormalStiffness << "\n";
    rLogger << "    Residual bond stress           : " << mResidualBondStress << "\n";
    rLogger << "    Maximum bond stress            : " << mMaxBondStress << "\n";
    rLogger << "    Alpha                          : " << mAlpha << "\n";
    rLogger << "    Slip at maximum bond stress    : " << mSlipAtMaxBondStress << "\n";
    rLogger << "    Slip at residual bond stress   : " << mSlipAtResidualBondStress << "\n";
}

// check parameters
void NuTo::FibreMatrixBondStressSlip::CheckParameters() const
{
    assert(mNormalStiffness > 0.0 and "Normal stiffnes is <= 0 or not initialized properly");
    assert(mResidualBondStress >= 0.0 and "Residual bond stress is <= 0 or not initialized properly");
    assert(mMaxBondStress > 0.0 and "Max bond stress is <= 0 or not initialized properly");
    assert(mSlipAtMaxBondStress > 0.0 and "Slip at max bond stress is <= 0 or not initialized properly");
    assert(mSlipAtResidualBondStress > 0.0 and "Slip at residual bond stress is <= 0 or not initialized properly");
}

namespace NuTo // template specialization in same namespace as definition
{
template<>
NuTo::eError NuTo::FibreMatrixBondStressSlip::Evaluate<1>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput,
    Data& staticData)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "IMPLEMENT ME!!!");
    return eError::SUCCESSFUL;
}

template<>
NuTo::eError NuTo::FibreMatrixBondStressSlip::Evaluate<2>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput,
    Data& rStaticData)
{
    const auto& slip = *rConstitutiveInput.at(Constitutive::eInput::INTERFACE_SLIP);
    const double slipLastConverged = GetCurrentStaticData(rStaticData, rConstitutiveInput);

    const double slipMaxHistory = std::max(slipLastConverged, slip(0, 0));

    bool performUpdateAtEnd = false;

    /////////////////////////////////////////////////
    //         LOOP OVER OUTPUT REQUESTS           //
    /////////////////////////////////////////////////

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
            case NuTo::Constitutive::eOutput::INTERFACE_CONSTITUTIVE_MATRIX:
            {
                ConstitutiveMatrixXd& constitutiveMatrix = dynamic_cast<ConstitutiveMatrixXd&>(*itOutput.second);
                constitutiveMatrix.setZero(mGlobalDimension, mGlobalDimension);

                // the first component of the constitutive matrix is the derivative of the interface stress with respect to the interface slip
                if (std::abs(slip(0, 0)) <= mSlipAtMaxBondStress)
                {
                    constitutiveMatrix(0, 0) =
                        mAlpha * mMaxBondStress * std::pow(slip(0, 0) / mSlipAtMaxBondStress, mAlpha - 1)
                            / mSlipAtMaxBondStress;

                }
                else if (std::abs(slip(0, 0)) > mSlipAtMaxBondStress
                    and std::abs(slip(0, 0)) <= mSlipAtResidualBondStress)
                {
                    constitutiveMatrix(0, 0) =
                        (mMaxBondStress - mResidualBondStress) / (mSlipAtMaxBondStress - mSlipAtResidualBondStress);

                }
                else
                {
                    constitutiveMatrix(0, 0) = 0.0;
                }

                // the bond stress-slip relationship needs to be adjusted if unloading occurs
                if (std::abs(slip(0, 0)) < slipMaxHistory)
                {
                    // first the bond stress-slipMaxHistory relationship is adjusted
                    double newMaxBondStress = 0.0;
                    if (std::abs(slipMaxHistory) <= mSlipAtMaxBondStress)
                    {
                        newMaxBondStress = mMaxBondStress * std::pow(slipMaxHistory / mSlipAtMaxBondStress, mAlpha);

                    }
                    else if (slipMaxHistory > mSlipAtMaxBondStress and slipMaxHistory <= mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = mMaxBondStress
                            - (mMaxBondStress - mResidualBondStress) * (slipMaxHistory - mSlipAtMaxBondStress)
                                / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

                    }
                    else if (slipMaxHistory < -mSlipAtMaxBondStress and slipMaxHistory >= -mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = -mMaxBondStress
                            - (-mMaxBondStress + mResidualBondStress) * (slipMaxHistory + mSlipAtMaxBondStress)
                                / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

                    }
                    else if (slipMaxHistory > mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = mResidualBondStress;

                    }
                    else if (slipMaxHistory < -mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = mResidualBondStress;
                    }
                    else
                    {
                        throw MechanicsException(std::string(__PRETTY_FUNCTION__)
                                                     + ":\t Check if clause. This branch should never be executed!");
                    }

                    // second the adjusted bond stress-slip relationship is evaluated
                    if (slip(0, 0) >= -slipMaxHistory and slip(0, 0) <= slipMaxHistory)
                    {
                        constitutiveMatrix(0, 0) =
                            mAlpha * newMaxBondStress * std::pow(slip(0, 0) / slipMaxHistory, mAlpha - 1)
                                / slipMaxHistory;

                    }
                    else if (slip(0, 0) > slipMaxHistory and slip(0, 0) <= mSlipAtResidualBondStress)
                    {
                        constitutiveMatrix(0, 0) =
                            (newMaxBondStress - mResidualBondStress) / (slipMaxHistory - mSlipAtResidualBondStress);

                    }
                    else if (slip(0, 0) < -slipMaxHistory and slip(0, 0) >= -mSlipAtResidualBondStress)
                    {
                        constitutiveMatrix(0, 0) =
                            (newMaxBondStress - mResidualBondStress) / (slipMaxHistory - mSlipAtResidualBondStress);

                    }
                    else
                    {
                        constitutiveMatrix(0, 0) = 0.0;
                    }
                }

                // the normal component of the interface slip is equal to the penalty stiffness to avoid penetration
                for (unsigned int iDim = 1; iDim < mGlobalDimension; ++iDim)
                    constitutiveMatrix(iDim, iDim) = mNormalStiffness;

                break;
            }

            case NuTo::Constitutive::eOutput::BOND_STRESS:
            {

                ConstitutiveMatrixXd& bondStress = dynamic_cast<ConstitutiveMatrixXd&>(*itOutput.second);
                bondStress.setZero(mGlobalDimension, 1);

                // the interface stresses correspond to a modified Bertero-EligeHausen-Popov bond stress-slip model
                if (std::abs(slip(0, 0)) <= mSlipAtMaxBondStress)
                {
                    bondStress(0, 0) = mMaxBondStress * std::pow(slip(0, 0) / mSlipAtMaxBondStress, mAlpha);

                }
                else if (slip(0, 0) > mSlipAtMaxBondStress and slip(0, 0) <= mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = mMaxBondStress
                        - (mMaxBondStress - mResidualBondStress) * (slip(0, 0) - mSlipAtMaxBondStress)
                            / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

                }
                else if (slip(0, 0) < -mSlipAtMaxBondStress and slip(0, 0) >= -mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = -mMaxBondStress
                        - (-mMaxBondStress + mResidualBondStress) * (slip(0, 0) + mSlipAtMaxBondStress)
                            / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

                }
                else if (slip(0, 0) > mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = mResidualBondStress;

                }
                else if (slip(0, 0) < -mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = -mResidualBondStress;

                }
                else
                {
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__)
                                                 + ":\t Check if clause. This branch should never be executed!");
                }

                // the bond stress-slip relationship needs to be adjusted if unloading occurs
                if (std::abs(slip(0, 0)) < slipMaxHistory)
                {
                    double newMaxBondStress = 0.0;
                    // first the bond stress-slip relationship is adjusted
                    if (slipMaxHistory >= -mSlipAtMaxBondStress and slipMaxHistory <= mSlipAtMaxBondStress)
                    {
                        newMaxBondStress = mMaxBondStress * std::pow(slipMaxHistory / mSlipAtMaxBondStress, mAlpha);

                    }
                    else if (slipMaxHistory > mSlipAtMaxBondStress and slipMaxHistory <= mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = mMaxBondStress
                            - (mMaxBondStress - mResidualBondStress) * (slipMaxHistory - mSlipAtMaxBondStress)
                                / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

                    }
                    else if (slipMaxHistory < -mSlipAtMaxBondStress and slipMaxHistory >= -mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = -mMaxBondStress
                            - (-mMaxBondStress + mResidualBondStress) * (slipMaxHistory + mSlipAtMaxBondStress)
                                / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

                    }
                    else if (slipMaxHistory > mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = mResidualBondStress;

                    }
                    else if (slipMaxHistory < -mSlipAtResidualBondStress)
                    {
                        newMaxBondStress = -mResidualBondStress;
                    }
                    else
                    {
                        throw MechanicsException(std::string(__PRETTY_FUNCTION__)
                                                     + ":\t Check if clause. This branch should never be executed!");
                    }

                    // second the adjusted bond stress-slip relationship is evaluated
                    if (slip(0, 0) >= -slipMaxHistory and slip(0, 0) <= slipMaxHistory)
                    {
                        bondStress(0, 0) = newMaxBondStress * std::pow(slip(0, 0) / slipMaxHistory, mAlpha);

                    }
                    else if (slip(0, 0) > slipMaxHistory and slip(0, 0) <= mSlipAtResidualBondStress)
                    {
                        bondStress(0, 0) = newMaxBondStress
                            - (newMaxBondStress - mResidualBondStress) * (slip(0, 0) - slipMaxHistory)
                                / (mSlipAtResidualBondStress - slipMaxHistory);

                    }
                    else if (slip(0, 0) < -slipMaxHistory and slip(0, 0) >= -mSlipAtResidualBondStress)
                    {
                        bondStress(0, 0) = -newMaxBondStress
                            - (-newMaxBondStress + mResidualBondStress) * (slip(0, 0) + slipMaxHistory)
                                / (-mSlipAtResidualBondStress + slipMaxHistory);

                    }
                    else if (slip(0, 0) > mSlipAtResidualBondStress)
                    {
                        bondStress(0, 0) = mResidualBondStress;

                    }
                    else if (slip(0, 0) < -mSlipAtResidualBondStress)
                    {
                        bondStress(0, 0) = -mResidualBondStress;
                    }
                    else
                    {
                        throw MechanicsException(std::string(__PRETTY_FUNCTION__)
                                                     + ":\t Check if clause. This branch should never be executed!");
                    }

                }

                // the normal component of the interface slip is equal to the penalty stiffness to avoid penetration
                for (unsigned int iDim = 1; iDim < mGlobalDimension; ++iDim)
                    bondStress(iDim, 0) = mNormalStiffness * slip(iDim, 0);

            }
                break;
            case NuTo::Constitutive::eOutput::SLIP:
                break;
            case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            {
                performUpdateAtEnd = true;
            }
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Output object "
                    + NuTo::Constitutive::OutputToString(itOutput.first) + " could not be calculated, check the "
                    "allocated material law and the section behavior.");
        }
        if(itOutput.second!=nullptr)
            itOutput.second->SetIsCalculated(true);
    }

    //update history variables
    if (performUpdateAtEnd)
        rStaticData.SetData(slipLastConverged);

    return eError::SUCCESSFUL;
}

template<>
NuTo::eError NuTo::FibreMatrixBondStressSlip::Evaluate<3>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput,
    Data& staticData)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "IMPLEMENT ME!!!");
}

}

NuTo::ConstitutiveInputMap NuTo::FibreMatrixBondStressSlip::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;
    constitutiveInputMap[Constitutive::eInput::INTERFACE_SLIP];
    return constitutiveInputMap;
}

double NuTo::FibreMatrixBondStressSlip::GetCurrentStaticData(Data& rStaticData,
        const ConstitutiveInputMap& rConstitutiveInput) const
{
    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::eInput::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "You need to specify the way the static data should be calculated (input list).");

    const auto& calculateStaticData = *static_cast<const ConstitutiveCalculateStaticData*>(itCalculateStaticData->second.get());

    switch (calculateStaticData.GetCalculateStaticData())
    {
        case eCalculateStaticData::USE_PREVIOUS:
        {
            return rStaticData.GetData();
        }

        case eCalculateStaticData::EULER_BACKWARD:
        {
            const auto& slip = *rConstitutiveInput.at(Constitutive::eInput::INTERFACE_SLIP);

            int index = calculateStaticData.GetIndexOfPreviousStaticData();
            double oldSlip = rStaticData.GetData(index);
            
            return std::max(slip(0,0), oldSlip);
        }

        case eCalculateStaticData::EULER_FORWARD:
        {
            auto itTimeStep = rConstitutiveInput.find(Constitutive::eInput::TIME_STEP);

            if (itTimeStep == rConstitutiveInput.end())
                throw MechanicsException(__PRETTY_FUNCTION__, "TimeStep input needed for EULER_FORWARD.");
            const auto& timeStep = *itTimeStep->second;

            assert(rStaticData.GetNumData() >= 2);

            return ConstitutiveCalculateStaticData::EulerForward(
                    rStaticData.GetData(1), rStaticData.GetData(2), timeStep);
        }

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Cannot calculate the static data in the requested way.");
    }
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::FibreMatrixBondStressSlip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::FibreMatrixBondStressSlip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::FibreMatrixBondStressSlip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::FibreMatrixBondStressSlip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::FibreMatrixBondStressSlip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::FibreMatrixBondStressSlip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::FibreMatrixBondStressSlip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize FibreMatrixBondStressSlip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
    & BOOST_SERIALIZATION_NVP(mMaxBondStress)
    & BOOST_SERIALIZATION_NVP(mResidualBondStress)
    & BOOST_SERIALIZATION_NVP(mSlipAtMaxBondStress)
    & BOOST_SERIALIZATION_NVP(mSlipAtResidualBondStress)
    & BOOST_SERIALIZATION_NVP(mNormalStiffness)
    & BOOST_SERIALIZATION_NVP(mAlpha);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize FibreMatrixBondStressSlip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::FibreMatrixBondStressSlip)
#endif // ENABLE_SERIALIZATION


bool NuTo::FibreMatrixBondStressSlip::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    return rTimeDerivative<1 && rDofRow == Node::eDof::DISPLACEMENTS && rDofCol ==Node::eDof::DISPLACEMENTS;
}
