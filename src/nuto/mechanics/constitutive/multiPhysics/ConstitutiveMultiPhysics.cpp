#include "ConstitutiveMultiPhysics.h"

#include "ConstitutiveStaticDataMultiPhysics.h"
#include <nuto/mechanics/elements/ElementBase.h>

#include <set>


//! @brief ... default constructor
NuTo::ConstitutiveMultiPhysics::ConstitutiveMultiPhysics()
    : ConstitutiveBase()
{}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase *NuTo::ConstitutiveMultiPhysics::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const NuTo::ElementBase *rElement) const
{
    return new ConstitutiveStaticDataMultiPhysics;
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase *NuTo::ConstitutiveMultiPhysics::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const NuTo::ElementBase *rElement) const
{
    return new ConstitutiveStaticDataMultiPhysics;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::ConstitutiveMultiPhysics::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    for (unsigned int i = 0; i < mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->CheckElementCompatibility(rElementType))
        {
            return true;
        }
    }
    return false;
}

//! @brief ... check parameters of the constitutive relationship
//! if one check fails, an exception is thrwon
void NuTo::ConstitutiveMultiPhysics::CheckParameters() const
{
    /*
    for (unsigned int i = 0; i < mConstitutiveLaws.size(); i++)
    {
        mConstitutiveLaws[i]->CheckParameters();
    }*/

    // PROTECTED! ---> Maybe friend class if needed
}

//! @brief ... evaluate the constitutive relation of every attached constitutive law in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::ConstitutiveMultiPhysics::Evaluate1D(NuTo::ElementBase *rElement,
                                                               int rIp,
                                                               const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase *> &rConstitutiveInput,
                                                               std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase *> &rConstitutiveOutput)
{
    for (unsigned int i = 0; i < mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->CheckElementCompatibility(rElement->GetEnumType()))
        {
            mConstitutiveLaws[i]->Evaluate1D(rElement,
                                             rIp,
                                             rConstitutiveInput,
                                             rConstitutiveOutput);
        }
    }
    return Error::SUCCESSFUL;
}


//! @brief ... evaluate the constitutive relation of every attached constitutive law in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::ConstitutiveMultiPhysics::Evaluate2D(NuTo::ElementBase *rElement,
                                                               int rIp,
                                                               const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase *> &rConstitutiveInput,
                                                               std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase *> &rConstitutiveOutput)
{
    std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase *> LawSpecificOutputs;
    std::set<NuTo::Constitutive::Output::eOutput> UnhandledOutputs;
    for (auto Output_it : rConstitutiveOutput)
    {
        UnhandledOutputs.insert(Output_it.first);
    }

    for (unsigned int i = 0; i < mConstitutiveLaws.size(); i++)
    {
         if(mConstitutiveLaws[i]->CheckElementCompatibility(rElement->GetEnumType()))
        {

            LawSpecificOutputs.clear();
            for (auto Output_it : rConstitutiveOutput)
            {
                if(mConstitutiveLaws[i]->CheckOutputTypeCompatibility(Output_it.first))
                {
                    LawSpecificOutputs.insert(Output_it);
                    auto UnhandledOutput_it = UnhandledOutputs.find(Output_it.first);
                    if (UnhandledOutput_it != UnhandledOutputs.end())
                    {
                        UnhandledOutputs.erase(UnhandledOutput_it);
                    }
                }
            }

            mConstitutiveLaws[i]->Evaluate2D(rElement,
                                             rIp,
                                             rConstitutiveInput,
                                             LawSpecificOutputs);
        }
    }
    if(!UnhandledOutputs.empty())
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMultiPhysics::Evaluate2D] There are constitutive outputs that aren't evaluated by any of the attached constitutive laws.");
    }

    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation of every attached constitutive law relation in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::ConstitutiveMultiPhysics::Evaluate3D(NuTo::ElementBase *rElement,
                                                               int rIp,
                                                               const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase *> &rConstitutiveInput,
                                                               std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase *> &rConstitutiveOutput)
{
    for (unsigned int i = 0; i < mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->CheckElementCompatibility(rElement->GetEnumType()))
        {
            mConstitutiveLaws[i]->Evaluate3D(rElement,
                                             rIp,
                                             rConstitutiveInput,
                                             rConstitutiveOutput);
        }
    }
    return Error::SUCCESSFUL;
}


//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::ConstitutiveMultiPhysics::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    bool ParameterFound = false;

    double ParameterReturn;
    double TempParam;
    for(unsigned int i=0; i<mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->CheckHaveParameter(rIdentifier))
        {
            if (ParameterFound)
            {
                TempParam = mConstitutiveLaws[i]->GetParameterDouble(rIdentifier);
                if (TempParam!=ParameterReturn)
                {
                    throw NuTo::MechanicsException("[NuTo::ConstitutiveMultiPhysics:GetParameterBool] A requested parameter is not consistent in all of the attached constitutive laws");
                }
            }
            else
            {
                ParameterReturn = mConstitutiveLaws[i]->GetParameterDouble(rIdentifier);
                ParameterFound = true;
            }
        }
    }
    if(!ParameterFound)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMultiPhysics::GetParameterBool] Tried to get a parameter that is not used by any of the attached constitutive laws.");
    }
    return ParameterReturn;
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable

void NuTo::ConstitutiveMultiPhysics::SetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    bool ParameterFound = false;
    for(unsigned int i=0; i<mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->CheckHaveParameter(rIdentifier))
        {
            mConstitutiveLaws[i]->SetParameterBool(rIdentifier,rValue);
            ParameterFound = true;
        }
    }
    if(!ParameterFound)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMultiPhysics::SetParameterBool] Tried to set a parameter that is not used by any of the attached constitutive laws.");
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveMultiPhysics::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    bool ParameterFound = false;
    for(unsigned int i=0; i<mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->CheckHaveParameter(rIdentifier))
        {
            mConstitutiveLaws[i]->SetParameterDouble(rIdentifier,rValue);
            ParameterFound = true;
        }
    }
    if(!ParameterFound)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMultiPhysics::SetParameterDouble] Tried to set a parameter that is not used by any of the attached constitutive laws.");
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::ConstitutiveMultiPhysics::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    bool ParameterFound = false;
    for(unsigned int i=0; i<mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->CheckHaveParameter(rIdentifier))
        {
            mConstitutiveLaws[i]->SetParameterFullVectorDouble(rIdentifier,rValue);
            ParameterFound = true;
        }
    }
    if(!ParameterFound)
    {
        throw NuTo::MechanicsException("[NuTo::ConstitutiveMultiPhysics::SetParameterFullVectorDouble] Tried to set a parameter that is not used by any of the attached constitutive laws.");
    }
}



//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::ConstitutiveMultiPhysics::GetType() const
{
    return Constitutive::eConstitutiveType::MULTI_PHYSICS;
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::ConstitutiveMultiPhysics::HaveTmpStaticData() const
{
    for (unsigned int i = 0; i < mConstitutiveLaws.size(); i++)
    {
        if(mConstitutiveLaws[i]->HaveTmpStaticData())
        {
            return true;
        }
    }
    return false;
}

//! @brief ... adds a constitutive law to a multi physics model
//! @param ... additional constitutive law
void NuTo::ConstitutiveMultiPhysics::MultiPhysicsAddConstitutiveLaw(NuTo::ConstitutiveBase *rConstitutiveLaw)
{
    mConstitutiveLaws.push_back(rConstitutiveLaw);
}
