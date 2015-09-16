#include "nuto/mechanics/constitutive/mechanics/InterfaceGoodman.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/mechanics/InterfaceSlip.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

NuTo::InterfaceGoodman::InterfaceGoodman() :
        ConstitutiveBase(), mNormalStiffness(0), mTangentialStiffness(0)
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

    // calculate engineering strain
    if (rConstitutiveInput.find(NuTo::Constitutive::Input::INTERFACE_SLIP) == rConstitutiveInput.end())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Interface slip needed to evaluate.");
    const InterfaceSlip& interfaceSlip(rConstitutiveInput.find(NuTo::Constitutive::Input::INTERFACE_SLIP)->second->GetInterfaceSlip());
    Eigen::VectorXd interfaceSlipVector = interfaceSlip.GetInterfaceSlipVector();

    CheckParameters();

    for (std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>::iterator itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch (itOutput->first)
        {
        case NuTo::Constitutive::Output::INTERFACE_CONSTITUTIVE_MATRIX:
        {
            ConstitutiveTangentLocal<2, 2>& constitutiveMatrix(itOutput->second->AsConstitutiveTangentLocal_2x2());

            constitutiveMatrix.setZero();
            constitutiveMatrix(0, 0) = mTangentialStiffness;
            constitutiveMatrix(1, 1) = mNormalStiffness;
        }
            break;
        case NuTo::Constitutive::Output::INTERFACE_STRESSES:
        {
            ConstitutiveTangentLocal<2, 1>& interfaceStresses(itOutput->second->AsConstitutiveTangentLocal_2x1());

            interfaceStresses.setZero();
            interfaceStresses(0, 0) = mTangentialStiffness * interfaceSlipVector(0,0);
            interfaceStresses(1, 0) = mNormalStiffness * interfaceSlipVector(1,0);
        }
            break;
        default:
            throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t output object)" + NuTo::Constitutive::OutputToString(itOutput->first) + std::string(" culd not be calculated, check the allocated material law and the section behavior."));
        }
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
    case Constitutive::eConstitutiveParameter::TANGENTIAL_STIFFNESS:
        return this->mTangentialStiffness;
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
    case Constitutive::eConstitutiveParameter::TANGENTIAL_STIFFNESS:
    {
        mTangentialStiffness = rValue;
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

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::InterfaceGoodman::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Normal stiffness               : " << this->mNormalStiffness << "\n";
    rLogger << "    Tangential stiffness           : " << this->mTangentialStiffness << "\n";
}

// check parameters
void NuTo::InterfaceGoodman::CheckParameters() const
{
    assert(mNormalStiffness > 0.0);
    assert(mTangentialStiffness > 0.0);
}

