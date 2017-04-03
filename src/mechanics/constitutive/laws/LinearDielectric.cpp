#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/laws/LinearDielectric.h"
#include "base/Logger.h"
#include "mechanics/MechanicsException.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

LinearDielectric::LinearDielectric() : ConstitutiveBase()
{
    mEps = 1.0;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
template <class Archive> void LinearDielectric::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LinearDielectric" << std::endl;
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase) & BOOST_SERIALIZATION_NVP(mK) &
            BOOST_SERIALIZATION_NVP(mCt) & BOOST_SERIALIZATION_NVP(mRho);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LinearDielectric" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(LinearDielectric)
#endif // ENABLE_SERIALIZATION

ConstitutiveInputMap LinearDielectric::GetConstitutiveInputs(
        const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType&) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        case Constitutive::eOutput::ELECTRIC_FIELD:
        {
            constitutiveInputMap[Constitutive::eInput::ELECTRIC_FIELD];
            break;
        }
        case Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        {
            // no input needed for the linear case
            // the result is just a constant
            break;
        }
        default:
            continue;
        }
    }

    return constitutiveInputMap;
}

bool LinearDielectric::CheckDofCombinationComputable(Node::eDof dofRow, Node::eDof dofCol, int timeDerivative) const
{
    if ((timeDerivative == 0) && (dofRow == Node::eDof::ELECTRICPOTENTIAL and dofCol == Node::eDof::ELECTRICPOTENTIAL))
    {
            return true;
    }
    else return false;
}

template <int TDim>
void LinearDielectric::Evaluate(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    auto eye = Eigen::MatrixXd::Identity(TDim, TDim);

    InputData<TDim> inputData;
    for (auto& itInput : rConstitutiveInput)
    {
        switch (itInput.first)
        {
        case Constitutive::eInput::ELECTRIC_FIELD:
            inputData.mElectricField = static_cast<ConstitutiveVector<TDim>*>(itInput.second.get())->AsVector();
            break;
        default:
            continue;
        }
    }

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::ELECTRIC_FIELD:
        {
            Eigen::Matrix<double, TDim, 1>& electricField = *static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get());
            electricField = inputData.mElectricField;
            break;
        }

        case Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        {
            Eigen::Matrix<double, TDim, 1>& electricDisplacement = *static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get());
            electricDisplacement = mEps * eye * inputData.mElectricField;
            break;
        }
        case Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        {
            Eigen::Matrix<double, TDim, TDim>& permittivity =
                    *static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second.get());
            permittivity = mEps * eye;
            break;
        }

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

bool LinearDielectric::CheckHaveParameter(Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::RELATIVE_ELECTRIC_PERMITTIVITY:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}

double LinearDielectric::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::RELATIVE_ELECTRIC_PERMITTIVITY:
        return this->mEps;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void LinearDielectric::SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
//    ConstitutiveBase::CheckParameterDouble(rIdentifier, rValue);
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::RELATIVE_ELECTRIC_PERMITTIVITY:
        this->mEps = rValue;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

bool LinearDielectric::CheckOutputTypeCompatibility(Constitutive::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}

Constitutive::eConstitutiveType LinearDielectric::GetType() const
{
    return Constitutive::eConstitutiveType::LINEAR_DIELECTRIC;
}

bool LinearDielectric::CheckElementCompatibility(Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case Element::eElementType::CONTINUUMELEMENT:
        return true;
    default:
        return false;
    }
}

void LinearDielectric::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Relative electric permittivity          : " << this->mEps << "\n";
}

void LinearDielectric::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::RELATIVE_ELECTRIC_PERMITTIVITY, mEps);
}

template void LinearDielectric::Evaluate<1>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void LinearDielectric::Evaluate<2>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void LinearDielectric::Evaluate<3>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
