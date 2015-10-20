#include "nuto/mechanics/constitutive/mechanics/InterfaceGoodman.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataBondStressSlip.h"
#include "nuto/mechanics/constitutive/mechanics/InterfaceSlip.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

NuTo::InterfaceGoodman::InterfaceGoodman() :
        ConstitutiveBase(), mMaxBondStress(0), mResidualBondStress(0), mSlipAtMaxBondStress(0), mSlipAtResidualBondStress(0), mAlpha(0), mNormalStiffness(0)
{
}

//! @brief ... evaluate the constitutive relation in 1D
NuTo::Error::eError NuTo::InterfaceGoodman::Evaluate1D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Not implemented.");
}

//! @brief ... evaluate the constitutive relation in 2D
NuTo::Error::eError NuTo::InterfaceGoodman::Evaluate2D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{

    CheckParameters();

    assert(mAlpha==1.0 and "Not implemented for arbitrary alpha yet");

    // interface slip
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::INTERFACE_SLIP) == rConstitutiveInput.end())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Interface slip needed to evaluate.");

    const InterfaceSlip& interfaceSlip(rConstitutiveInput.find(NuTo::Constitutive::Input::INTERFACE_SLIP)->second->GetInterfaceSlip());

    Eigen::VectorXd interfaceSlipVector = interfaceSlip.GetInterfaceSlipVector();

    //Get previous ip_data
    ConstitutiveStaticDataBondStressSlip *oldSlip = (rElement->GetStaticData(rIp))->AsBondStressSlip();

    const double slipLastConverged = oldSlip->GetSlip();
    const double slip = std::max(slipLastConverged, std::abs(interfaceSlipVector.at(0, 0)));

    //set this to true, if update is in the map, perform the update after all other outputs have been calculated
    bool performUpdateAtEnd(false);

    /////////////////////////////////////////////////
    //         LOOP OVER OUTPUT REQUESTS           //
    /////////////////////////////////////////////////

    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch (itOutput->first)
        {
        case NuTo::Constitutive::Output::INTERFACE_CONSTITUTIVE_MATRIX:
        {
            ConstitutiveTangentLocal<2, 2>& constitutiveMatrix(itOutput->second->AsConstitutiveTangentLocal_2x2());

            constitutiveMatrix.setZero();

            // the first component of the constitutive matrix is the derivative of the interface stress with respect to the interface slip
            if (std::abs(interfaceSlipVector.at(0, 0)) <= mSlipAtMaxBondStress)
            {
                constitutiveMatrix(0, 0) = mAlpha * mMaxBondStress * std::pow(interfaceSlipVector.at(0, 0) / mSlipAtMaxBondStress, mAlpha - 1) / mSlipAtMaxBondStress;

            } else if (interfaceSlipVector.at(0, 0) > mSlipAtMaxBondStress and interfaceSlipVector.at(0, 0) <= mSlipAtResidualBondStress)
            {
                constitutiveMatrix(0, 0) = (mMaxBondStress - mResidualBondStress) / (mSlipAtMaxBondStress - mSlipAtResidualBondStress);

            } else if (interfaceSlipVector.at(0, 0) < -mSlipAtMaxBondStress and interfaceSlipVector.at(0, 0) >= -mSlipAtResidualBondStress)
            {
                constitutiveMatrix(0, 0) = (mMaxBondStress - mResidualBondStress) / (mSlipAtMaxBondStress - mSlipAtResidualBondStress);

            } else
            {
                constitutiveMatrix(0, 0) = 0.0;
            }


            // the normal component of the interface slip is equal to the penalty stiffness to avoid penetration
            constitutiveMatrix(1, 1) = mNormalStiffness;


            // the bond stress-slip relationship needs to be adjusted if unloading occurs
            if (std::abs(interfaceSlipVector.at(0, 0)) < slip)
            {
                // first the bond stress-slip relationship is adjusted
                double newMaxBondStress = 0.0;
                if (std::abs(slip) <= mSlipAtMaxBondStress)
                {
                    newMaxBondStress = mMaxBondStress * std::pow(slip / mSlipAtMaxBondStress, mAlpha);

                } else if (slip > mSlipAtMaxBondStress and slip <= mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mMaxBondStress - (mMaxBondStress - mResidualBondStress) * (slip - mSlipAtMaxBondStress) / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

                } else if (slip < -mSlipAtMaxBondStress and slip >= -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = -mMaxBondStress - (-mMaxBondStress + mResidualBondStress) * (slip + mSlipAtMaxBondStress) / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

                } else if (slip > mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mResidualBondStress;

                } else if (slip < -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mResidualBondStress;
                } else
                {
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
                }

                // second the adjusted bond stress-slip relationship is evaluated
                if (interfaceSlipVector(0, 0) >= -slip and interfaceSlipVector(0, 0) <= slip)
                {
                    constitutiveMatrix(0, 0) = mAlpha * newMaxBondStress * std::pow(interfaceSlipVector.at(0, 0) / slip, mAlpha - 1) / slip;

                } else if (interfaceSlipVector.at(0, 0) > slip and interfaceSlipVector.at(0, 0) <= mSlipAtResidualBondStress)
                {
                    constitutiveMatrix(0, 0) = (newMaxBondStress - mResidualBondStress) / (slip - mSlipAtResidualBondStress);

                } else if (interfaceSlipVector.at(0, 0) < -slip and interfaceSlipVector.at(0, 0) >= -mSlipAtResidualBondStress)
                {
                    constitutiveMatrix(0, 0) = (newMaxBondStress - mResidualBondStress) / (slip - mSlipAtResidualBondStress);

                } else
                {
                    constitutiveMatrix(0, 0) = 0.0;
                }
            }


        }
            break;
        case NuTo::Constitutive::Output::INTERFACE_STRESSES:
        {
            ConstitutiveTangentLocal<2, 1>& interfaceStresses(itOutput->second->AsConstitutiveTangentLocal_2x1());

            interfaceStresses.setZero();

            // the interface stresses correspond to a modified Bertero-EligeHausen-Popov bond stress-slip model
            if (std::abs(interfaceSlipVector.at(0, 0)) <= mSlipAtMaxBondStress)
            {
                interfaceStresses(0, 0) = mMaxBondStress * std::pow(interfaceSlipVector.at(0, 0) / mSlipAtMaxBondStress, mAlpha);

            } else if (interfaceSlipVector.at(0, 0) > mSlipAtMaxBondStress and interfaceSlipVector.at(0, 0) <= mSlipAtResidualBondStress)
            {
                interfaceStresses(0, 0) = mMaxBondStress - (mMaxBondStress - mResidualBondStress) * (interfaceSlipVector(0, 0) - mSlipAtMaxBondStress) / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

            } else if (interfaceSlipVector.at(0, 0) < -mSlipAtMaxBondStress and interfaceSlipVector.at(0, 0) >= -mSlipAtResidualBondStress)
            {
                interfaceStresses(0, 0) = -mMaxBondStress - (-mMaxBondStress + mResidualBondStress) * (interfaceSlipVector(0, 0) + mSlipAtMaxBondStress) / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

            } else if (interfaceSlipVector.at(0, 0) > mSlipAtResidualBondStress or interfaceSlipVector.at(0, 0) < -mSlipAtResidualBondStress)
            {
                interfaceStresses(0, 0) = mResidualBondStress;

            } else
            {
                throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
            }



            // the bond stress-slip relationship needs to be adjusted if unloading occurs
            if (std::abs(interfaceSlipVector.at(0,0)) < slip)
            {
                double newMaxBondStress = 0.0;
                // first the bond stress-slip relationship is adjusted
                if (slip >= -mSlipAtMaxBondStress and slip <= mSlipAtMaxBondStress)
                {
                    newMaxBondStress = mMaxBondStress * std::pow(slip / mSlipAtMaxBondStress, mAlpha);

                } else if (slip > mSlipAtMaxBondStress and slip <= mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mMaxBondStress - (mMaxBondStress - mResidualBondStress) * (slip - mSlipAtMaxBondStress) / (mSlipAtResidualBondStress - mSlipAtMaxBondStress);

                } else if (slip < -mSlipAtMaxBondStress and slip >= -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = -mMaxBondStress - (-mMaxBondStress + mResidualBondStress) * (slip + mSlipAtMaxBondStress) / (-mSlipAtResidualBondStress + mSlipAtMaxBondStress);

                } else if (slip > mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mResidualBondStress;

                } else if (slip < -mSlipAtResidualBondStress)
                {
                    newMaxBondStress = mResidualBondStress;
                } else
                {
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
                }

                // second the adjusted bond stress-slip relationship is evaluated
                if (interfaceSlipVector(0, 0) >= -slip and interfaceSlipVector(0, 0) <= slip)
                {
                    interfaceStresses(0, 0) = newMaxBondStress * std::pow(interfaceSlipVector.at(0, 0) / slip, mAlpha);

                } else if (interfaceSlipVector(0, 0) > slip and interfaceSlipVector(0, 0) <= mSlipAtResidualBondStress)
                {
                    interfaceStresses(0, 0) = newMaxBondStress - (newMaxBondStress - mResidualBondStress) * (interfaceSlipVector(0, 0) - slip) / (mSlipAtResidualBondStress - slip);

                } else if (interfaceSlipVector(0, 0) < -slip and interfaceSlipVector(0, 0) >= -mSlipAtResidualBondStress)
                {
                    interfaceStresses(0, 0) = -newMaxBondStress - (-newMaxBondStress + mResidualBondStress) * (interfaceSlipVector(0, 0) + slip) / (-mSlipAtResidualBondStress + slip);

                } else if (interfaceSlipVector(0, 0) > mSlipAtResidualBondStress)
                {
                    interfaceStresses(0, 0) = mResidualBondStress;

                } else if (interfaceSlipVector(0, 0) < -mSlipAtResidualBondStress)
                {
                    interfaceStresses(0, 0) = mResidualBondStress;
                } else
                {
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Check if clause. This branch should never be executed!");
                }

            }



            interfaceStresses(1, 0) = mNormalStiffness * interfaceSlipVector(1, 0);

        }
            break;
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            break;
        default:
            throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t output object)" + NuTo::Constitutive::OutputToString(itOutput->first) + std::string(" culd not be calculated, check the allocated material law and the section behavior."));
        }
    }

    //update history variables
    if (performUpdateAtEnd)
    {
        oldSlip->SetSlip(slip);
    }
    return Error::SUCCESSFUL;

}

//! @brief ... evaluate the constitutive relation in 3D
NuTo::Error::eError NuTo::InterfaceGoodman::Evaluate3D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Not implemented.");
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
double NuTo::InterfaceGoodman::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
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
void NuTo::InterfaceGoodman::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
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
NuTo::Constitutive::eConstitutiveType NuTo::InterfaceGoodman::GetType() const
{
    return NuTo::Constitutive::INTERFACE_GOODMAN;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
bool NuTo::InterfaceGoodman::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::ELEMENT2DINTERFACE:
        return true;
    default:
        return false;
    }
}

NuTo::ConstitutiveStaticDataBase* NuTo::InterfaceGoodman::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataBondStressSlip();
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::InterfaceGoodman::Info(unsigned short rVerboseLevel, Logger& rLogger) const
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
void NuTo::InterfaceGoodman::CheckParameters() const
{
    assert(mNormalStiffness > 0.0 and "Normal stiffnes is <= 0 or not initialized properly");
    assert(mResidualBondStress >= 0.0 and "Residual bond stress is <= 0 or not initialized properly");
    assert(mMaxBondStress > 0.0 and "Max bond stress is <= 0 or not initialized properly");
    assert(mSlipAtMaxBondStress > 0.0 and "Slip at max bond stress is <= 0 or not initialized properly");
    assert(mSlipAtResidualBondStress > 0.0 and "Slip at residual bond stress is <= 0 or not initialized properly");
}

