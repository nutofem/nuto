// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/laws/ContactConstitutiveLaw.h"
#include "nuto/mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"

NuTo::ContactConstitutiveLaw::ContactConstitutiveLaw() :
        ConstitutiveBase()
{
}

NuTo::ContactConstitutiveLaw::ContactConstitutiveLaw(const std::function<double(double)> &rFunction, const std::function<double(double)> &rFunctionDerivative) :
        ConstitutiveBase(), mFunction(rFunction), mFunctionDerivative(rFunctionDerivative)
{
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::ContactConstitutiveLaw::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ContactConstitutiveLaw" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
    & BOOST_SERIALIZATION_NVP(mFunction)
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ContactConstitutiveLaw" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ContactConstitutiveLaw)
#endif // ENABLE_SERIALIZATION



NuTo::ConstitutiveInputMap NuTo::ContactConstitutiveLaw::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];

            if (rInterpolationType.IsConstitutiveInput(Node::eDof::TEMPERATURE))
                constitutiveInputMap[Constitutive::eInput::TEMPERATURE];

            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
            break;
        // no inputs needed for:
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            continue;
//            ProcessUnhandledOutput(__PRETTY_FUNCTION__,itOutput.first);
//            throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] output object " + Constitutive::OutputToString(itOutput.first) + " cannot be calculated by this constitutive law.");
        }
    }

    return constitutiveInputMap;
}

double NuTo::ContactConstitutiveLaw::GetContactForce(double rGap) const
{
    return mFunction(rGap);
}


double NuTo::ContactConstitutiveLaw::GetContactForceDerivative(double rGap) const
{
    return mFunctionDerivative(rGap);
}

bool NuTo::ContactConstitutiveLaw::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow,
                                                  NuTo::Node::eDof rDofCol,
                                                  int rTimeDerivative) const
{
    assert(rTimeDerivative == 0 or rTimeDerivative == 1 or rTimeDerivative == 2);
    switch (rTimeDerivative)
    {
    case 0:
    case 2:
    {
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            return true;
        default:
            return false;
        }

    }
        break;
    default:
        return false;
    }
}

bool NuTo::ContactConstitutiveLaw::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION:
        {
            return true;
        }
        case Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION:
        {
            return true;
        }
        default:
        {
            return false;
        }
    }
}

void NuTo::ContactConstitutiveLaw::SetParameterFunction(Constitutive::eConstitutiveParameter rIdentifier, const std::function<double(double)> &rFunction)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION:
        mFunction = rFunction;
        break;
    case Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION:
        mFunctionDerivative = rFunction;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
    SetParametersValid();
}

bool NuTo::ContactConstitutiveLaw::CheckOutputTypeCompatibility(NuTo::Constitutive::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    default:
    {
        return false;
    }
    }
}



NuTo::Constitutive::eConstitutiveType NuTo::ContactConstitutiveLaw::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::CONTACT_CONSTITUTIVE_LAW;
}



bool NuTo::ContactConstitutiveLaw::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
    case NuTo::Element::eElementType::CONTINUUMELEMENTIGA:
    case NuTo::Element::eElementType::CONTINUUMCONTACTELEMENT:
        return true;
    default:
        return false;
    }
}


void NuTo::ContactConstitutiveLaw::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    mFunction               : \n";
    rLogger << "    mFunctionDerivative     : \n";
}



void NuTo::ContactConstitutiveLaw::CheckParameters() const
{

}



