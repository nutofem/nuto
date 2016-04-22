#include "nuto/mechanics/constitutive/laws/FibreMatrixBondStressSlip.h"

//#include "nuto/mechanics/constitutive/mechanics/InterfaceSlip.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBondStressSlip.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrixXd.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/elements/IpDataStaticDataBase.h"

NuTo::FibreMatrixBondStressSlip::FibreMatrixBondStressSlip() :
        ConstitutiveBase(),
        mMaxBondStress(0),
        mResidualBondStress(0),
        mSlipAtMaxBondStress(0),
        mSlipAtResidualBondStress(0),
        mAlpha(0),
        mNormalStiffness(0)
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
    return NuTo::Constitutive::FIBRE_MATRIX_BOND_STRESS_SLIP;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
bool NuTo::FibreMatrixBondStressSlip::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::ELEMENT2DINTERFACE:
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

NuTo::Error::eError NuTo::FibreMatrixBondStressSlip::Evaluate1D(ElementBase* rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "IMPLEMENT ME!!!");
    return Error::SUCCESSFUL;
}

NuTo::Error::eError NuTo::FibreMatrixBondStressSlip::Evaluate2D(ElementBase* rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    const unsigned globalDimension = rElement->GetStructure()->GetDimension();
    const auto& slip = *rConstitutiveInput.at(Constitutive::Input::INTERFACE_SLIP);
    ConstitutiveStaticDataBondStressSlip currentStaticData = GetCurrentStaticData(*rElement, rIp, rConstitutiveInput);

    const double slipLastConverged = currentStaticData.GetSlip();
    const double slipMaxHistory = std::max(slipLastConverged, slip(0,0));

    bool performUpdateAtEnd = false;



    /////////////////////////////////////////////////
    //         LOOP OVER OUTPUT REQUESTS           //
    /////////////////////////////////////////////////

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::INTERFACE_CONSTITUTIVE_MATRIX:
        {
            ConstitutiveMatrixXd& constitutiveMatrix = dynamic_cast<ConstitutiveMatrixXd&>(*itOutput.second);
            constitutiveMatrix.setZero(globalDimension,globalDimension);

            // the first component of the constitutive matrix is the derivative of the interface stress with respect to the interface slip
            if (std::abs(slip(0,0)) <= mSlipAtMaxBondStress)
            {
                constitutiveMatrix(0, 0) = mAlpha * mMaxBondStress * std::pow(slip(0,0) / mSlipAtMaxBondStress, mAlpha - 1) / mSlipAtMaxBondStress;

            } else if (std::abs(slip(0,0)) > mSlipAtMaxBondStress and std::abs(slip(0,0)) <= mSlipAtResidualBondStress)
            {
                constitutiveMatrix(0, 0) = (mMaxBondStress - mResidualBondStress) / (mSlipAtMaxBondStress - mSlipAtResidualBondStress);

            } else
            {
                constitutiveMatrix(0, 0) = 0.0;
            }

            // the bond stress-slip relationship needs to be adjusted if unloading occurs
            if (std::abs(slip(0,0)) < slipMaxHistory)
            {
                // first the bond stress-slipMaxHistory relationship is adjusted
                double newMaxBondStress = 0.0;
                if (std::abs(slipMaxHistory) <= mSlipAtMaxBondStress)
                {
                    newMaxBondStress = mMaxBondStress * std::pow(slipMaxHistory / mSlipAtMaxBondStress, mAlpha);

                } else if (slipMaxHistory > mSlipAtMaxBondStress and slipMaxHistory <= mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mMaxBondStress - (mMaxBondStress - mResidualBondStress) * (slipMaxHistory - mSlipAtMaxBondStress) / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

                } else if (slipMaxHistory < -mSlipAtMaxBondStress and slipMaxHistory >= -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = -mMaxBondStress - (-mMaxBondStress + mResidualBondStress) * (slipMaxHistory + mSlipAtMaxBondStress) / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

                } else if (slipMaxHistory > mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mResidualBondStress;

                } else if (slipMaxHistory < -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mResidualBondStress;
                } else
                {
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
                }

                // second the adjusted bond stress-slip relationship is evaluated
                if (slip(0,0) >= -slipMaxHistory and slip(0,0) <= slipMaxHistory)
                {
                    constitutiveMatrix(0, 0) = mAlpha * newMaxBondStress * std::pow(slip(0,0) / slipMaxHistory, mAlpha - 1) / slipMaxHistory;

                } else if (slip(0,0) > slipMaxHistory and slip(0,0) <= mSlipAtResidualBondStress)
                {
                    constitutiveMatrix(0, 0) = (newMaxBondStress - mResidualBondStress) / (slipMaxHistory - mSlipAtResidualBondStress);

                } else if (slip(0,0) < -slipMaxHistory and slip(0,0) >= -mSlipAtResidualBondStress)
                {
                    constitutiveMatrix(0, 0) = (newMaxBondStress - mResidualBondStress) / (slipMaxHistory - mSlipAtResidualBondStress);

                } else
                {
                    constitutiveMatrix(0, 0) = 0.0;
                }
            }

            // the normal component of the interface slip is equal to the penalty stiffness to avoid penetration
            for (unsigned int iDim = 1; iDim < globalDimension; ++iDim)
                constitutiveMatrix(iDim, iDim) = mNormalStiffness;

            break;
        }

        case NuTo::Constitutive::Output::BOND_STRESS:
        {

            ConstitutiveMatrixXd& bondStress = dynamic_cast<ConstitutiveMatrixXd&>(*itOutput.second);
            bondStress.setZero(globalDimension,1);

            // the interface stresses correspond to a modified Bertero-EligeHausen-Popov bond stress-slip model
            if (std::abs(slip(0,0)) <= mSlipAtMaxBondStress)
            {
                bondStress(0, 0) = mMaxBondStress * std::pow(slip(0,0) / mSlipAtMaxBondStress, mAlpha);

            } else if (slip(0,0) > mSlipAtMaxBondStress and slip(0,0) <= mSlipAtResidualBondStress)
            {
                bondStress(0, 0) = mMaxBondStress - (mMaxBondStress - mResidualBondStress) * (slip(0,0) - mSlipAtMaxBondStress) / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

            } else if (slip(0,0) < -mSlipAtMaxBondStress and slip(0,0) >= -mSlipAtResidualBondStress)
            {
                bondStress(0, 0) = -mMaxBondStress - (-mMaxBondStress + mResidualBondStress) * (slip(0,0) + mSlipAtMaxBondStress) / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

            } else if (slip(0,0) > mSlipAtResidualBondStress)
            {
                bondStress(0, 0) = mResidualBondStress;

            } else if (slip(0,0) < -mSlipAtResidualBondStress)
            {
                bondStress(0, 0) = -mResidualBondStress;

            } else
            {
                throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
            }

            // the bond stress-slip relationship needs to be adjusted if unloading occurs
            if (std::abs(slip(0,0)) < slipMaxHistory)
            {
                double newMaxBondStress = 0.0;
                // first the bond stress-slip relationship is adjusted
                if (slipMaxHistory >= -mSlipAtMaxBondStress and slipMaxHistory <= mSlipAtMaxBondStress)
                {
                    newMaxBondStress = mMaxBondStress * std::pow(slipMaxHistory / mSlipAtMaxBondStress, mAlpha);

                } else if (slipMaxHistory > mSlipAtMaxBondStress and slipMaxHistory <= mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mMaxBondStress - (mMaxBondStress - mResidualBondStress) * (slipMaxHistory - mSlipAtMaxBondStress) / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

                } else if (slipMaxHistory < -mSlipAtMaxBondStress and slipMaxHistory >= -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = -mMaxBondStress - (-mMaxBondStress + mResidualBondStress) * (slipMaxHistory + mSlipAtMaxBondStress) / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

                } else if (slipMaxHistory > mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mResidualBondStress;

                } else if (slipMaxHistory < -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = -mResidualBondStress;
                } else
                {
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
                }

                // second the adjusted bond stress-slip relationship is evaluated
                if (slip(0,0)>= -slipMaxHistory and slip(0,0) <= slipMaxHistory)
                {
                    bondStress(0, 0) = newMaxBondStress * std::pow(slip(0,0) / slipMaxHistory, mAlpha);

                } else if (slip(0,0) > slipMaxHistory and slip(0,0) <= mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = newMaxBondStress - (newMaxBondStress - mResidualBondStress) * (slip(0,0) - slipMaxHistory) / (mSlipAtResidualBondStress - slipMaxHistory);

                } else if (slip(0,0) < -slipMaxHistory and slip(0,0) >= -mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = -newMaxBondStress - (-newMaxBondStress + mResidualBondStress) * (slip(0,0) + slipMaxHistory) / (-mSlipAtResidualBondStress + slipMaxHistory);

                } else if (slip(0,0) > mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = mResidualBondStress;

                } else if (slip(0,0) < -mSlipAtResidualBondStress)
                {
                    bondStress(0, 0) = -mResidualBondStress;
                } else
                {
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
                }

            }

            // the normal component of the interface slip is equal to the penalty stiffness to avoid penetration
            for (unsigned int iDim = 1; iDim < globalDimension; ++iDim)
                bondStress(iDim, 0) = mNormalStiffness * slip(iDim,0);


        }
            break;
        case NuTo::Constitutive::Output::SLIP:
            break;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            break;
        default:
            throw MechanicsException(std::string(__PRETTY_FUNCTION__, "output object ") + NuTo::Constitutive::OutputToString(itOutput.first) + std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    //update history variables
    if (performUpdateAtEnd)
        rElement->GetStaticData(rIp)->AsBondStressSlip()->SetSlip(currentStaticData.GetSlip());

    return Error::SUCCESSFUL;
}

NuTo::Error::eError NuTo::FibreMatrixBondStressSlip::Evaluate3D(ElementBase* rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "IMPLEMENT ME!!!");
    return Error::SUCCESSFUL;
}

NuTo::ConstitutiveStaticDataBase* NuTo::FibreMatrixBondStressSlip::AllocateStaticData1D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataBondStressSlip;
}

NuTo::ConstitutiveStaticDataBase* NuTo::FibreMatrixBondStressSlip::AllocateStaticData2D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataBondStressSlip;
}

NuTo::ConstitutiveStaticDataBase* NuTo::FibreMatrixBondStressSlip::AllocateStaticData3D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataBondStressSlip;
}

NuTo::ConstitutiveInputMap NuTo::FibreMatrixBondStressSlip::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;
    constitutiveInputMap[Constitutive::Input::INTERFACE_SLIP];
    return constitutiveInputMap;
}

NuTo::ConstitutiveStaticDataBondStressSlip NuTo::FibreMatrixBondStressSlip::GetCurrentStaticData(ElementBase& rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput) const
{
    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::Input::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "You need to specify the way the static data should be calculated (input list).");

    const auto& calculateStaticData = dynamic_cast<const ConstitutiveCalculateStaticData&>(*itCalculateStaticData->second);

    switch (calculateStaticData.GetCalculateStaticData())
    {
        case CalculateStaticData::USE_PREVIOUS:
        {
            return *(rElement.GetStaticData(rIp)->AsBondStressSlip());
        }

        case CalculateStaticData::EULER_BACKWARD:
        {
            const auto& slip = *rConstitutiveInput.at(Constitutive::Input::INTERFACE_SLIP);

            int index = calculateStaticData.GetIndexOfPreviousStaticData();
            const ConstitutiveStaticDataBondStressSlip& oldStaticData = *(rElement.GetStaticDataBase(rIp).GetStaticData(index)->AsBondStressSlip());

            ConstitutiveStaticDataBondStressSlip newStaticData;
            newStaticData.SetSlip(std::max(slip(0,0), oldStaticData.GetSlip()));

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

            ConstitutiveStaticDataBondStressSlip newStaticData;
            double newSlip = ConstitutiveCalculateStaticData::EulerForward(
                    staticData.GetStaticData(1)->AsBondStressSlip()->GetSlip(),
                    staticData.GetStaticData(2)->AsBondStressSlip()->GetSlip(),
                    timeStep);

//            std::cout << newKappa << std::endl;

            newStaticData.SetSlip(newSlip);
            return newStaticData;
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


bool NuTo::FibreMatrixBondStressSlip::CheckDofCombinationComputeable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    assert(rTimeDerivative>-1);
    if(rTimeDerivative<1 &&
       rDofRow == Node::DISPLACEMENTS &&
       rDofCol ==Node::DISPLACEMENTS)
    {
        return true;
    }
    return false;
}
