// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "base/Logger.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveInputBase.h"
#include "mechanics/constitutive/ConstitutiveOutputBase.h"
#include "mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "mechanics/constitutive/mechanics/Damage.h"
#include "mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "mechanics/constitutive/mechanics/LinearSpring.h"
#include "mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/sections/SectionBase.h"
#include "mechanics/sections/SectionEnum.h"

NuTo::LinearSpring::LinearSpring() :
        ConstitutiveBase(), mSpringStiffness(0)
{
    mSpringDirection.setZero();
}

//! @brief ... evaluate the constitutive relation in 1D
NuTo::Error::eError NuTo::LinearSpring::Evaluate1D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    return Error::SUCCESSFUL;
}

//! @brief ... evaluate the constitutive relation in 2D
NuTo::Error::eError NuTo::LinearSpring::Evaluate2D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw MechanicsException("[NuTo::LinearSpring::Evaluate2D] Not implemented.");
}

//! @brief ... evaluate the constitutive relation in 3D
NuTo::Error::eError NuTo::LinearSpring::Evaluate3D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    throw MechanicsException("[NuTo::LinearSpring::Evaluate3D] Not implemented.");
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::LinearSpring::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const
{
    return 0;
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::LinearSpring::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const
{
    return 0;
}

//! @brief ... allocate the correct static data
//! @return ... see brief explanation
NuTo::ConstitutiveStaticDataBase* NuTo::LinearSpring::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const
{
    return 0;
}

///////////////////////////////////////////////////////////////////////////

// parameters /////////////////////////////////////////////////////////////

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
double NuTo::LinearSpring::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    CheckParameters();

    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::SPRING_STIFFNESS:
        return this->mSpringStiffness;
    default:
    {
        throw MechanicsException("[NuTo::LinearSpring::GetParameterDouble] Constitutive law does not have the requested variable");
    }
    }
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::LinearSpring::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::SPRING_STIFFNESS:
    {
        mSpringStiffness = rValue;
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::LinearSpring::SetParameterDouble] Constitutive law does not have the requested variable");
    }
    }

}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::LinearSpring::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    CheckParameters();

    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::SPRING_DIRECTION:
    {
        return this->mSpringDirection;
    }
    default:
    {
        throw MechanicsException("[NuTo::LinearSpring::GetParameterFullVectorDouble] Constitutive law does not have the requested variable");
    }
    }
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::LinearSpring::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::SPRING_DIRECTION:
    {
        mSpringDirection = rValue;
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::LinearSpring::SetParameterFullVectorDouble] Constitutive law does not have the requested variable");
    }
    }

}

///////////////////////////////////////////////////////////////////////////

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::LinearSpring::GetType() const
{
    return NuTo::Constitutive::LINEAR_SPRING;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::LinearSpring::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::ELEMENT1DSPRING:
        return true;
    default:
        return false;
    }
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::LinearSpring::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Spring stiffness               : " << this->mSpringStiffness << "\n";
    rLogger << "    Spring direction               : " << this->mSpringDirection << "\n";
}

// check parameters
void NuTo::LinearSpring::CheckParameters() const
{
    assert(mSpringStiffness > 0.0);
    assert(mSpringDirection.norm() > 1e-8);
}

