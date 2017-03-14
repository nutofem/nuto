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
#include "mechanics/constitutive/laws/HeatConduction.h"
#include "base/Logger.h"
#include "mechanics/MechanicsException.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

HeatConduction::HeatConduction() : ConstitutiveBase()
{
    mK = 0.0;
    mCt = 0.0;
    mRho = 0.0;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
template <class Archive> void HeatConduction::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize HeatConduction" << std::endl;
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase) & BOOST_SERIALIZATION_NVP(mK) &
            BOOST_SERIALIZATION_NVP(mCt) & BOOST_SERIALIZATION_NVP(mRho);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize HeatConduction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(HeatConduction)
#endif // ENABLE_SERIALIZATION

ConstitutiveInputMap HeatConduction::GetConstitutiveInputs(
        const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType&) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::HEAT_FLUX:
        {
            constitutiveInputMap[Constitutive::eInput::TEMPERATURE_GRADIENT];
            break;
        }
        case Constitutive::eOutput::HEAT_CHANGE:
        {
            constitutiveInputMap[Constitutive::eInput::TEMPERATURE_CHANGE];
            break;
        }
        case Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        case Constitutive::eOutput::D_HEAT_D_TEMPERATURE:
        {
            // for nonlinear:
            // constitutiveInputMap[Constitutive::eInput::TEMPERATURE];
            break;
        }
        case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            continue;
        }
    }

    return constitutiveInputMap;
}

bool HeatConduction::CheckDofCombinationComputable(Node::eDof dofRow, Node::eDof dofCol, int timeDerivative) const
{
    if (timeDerivative == 2) return false;
    else if (dofRow == Node::eDof::TEMPERATURE and dofCol == Node::eDof::TEMPERATURE) return true;
    else return false;
}

template <int TDim>
void HeatConduction::Evaluate(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    auto eye = Eigen::MatrixXd::Identity(TDim, TDim);

    InputData<TDim> inputData;
    for (auto& itInput : rConstitutiveInput)
    {
        switch (itInput.first)
        {
        case Constitutive::eInput::TEMPERATURE_GRADIENT:
            inputData.mTemperatureGradient = static_cast<ConstitutiveVector<TDim>*>(itInput.second.get())->AsVector();
            break;
        case Constitutive::eInput::TEMPERATURE_CHANGE:
        {
            auto tempChange = *static_cast<ConstitutiveScalar*>(itInput.second.get());
            inputData.mTemperatureChange = tempChange[0];
            break;
        }
        case Constitutive::eInput::CALCULATE_STATIC_DATA:
        case Constitutive::eInput::TIME_STEP:
            break;
        default:
            continue;
        }
    }

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::HEAT_FLUX:
        {
            Eigen::Matrix<double, TDim, 1>& heatFlux = *static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get());
            heatFlux = -mK * eye * inputData.mTemperatureGradient;
            break;
        }
        case Constitutive::eOutput::HEAT_CHANGE:
        {
            Eigen::Matrix<double, 1, 1>& heatChange = *static_cast<ConstitutiveScalar*>(itOutput.second.get());
            heatChange(0, 0) = mCt * mRho * inputData.mTemperatureChange;
            break;
        }
        case Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        {
            Eigen::Matrix<double, TDim, TDim>& conductivity =
                    *static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second.get());
            conductivity = mK * eye;
            break;
        }
        case Constitutive::eOutput::D_HEAT_D_TEMPERATURE:
        {
            Eigen::Matrix<double, 1, 1>& tangent = *static_cast<ConstitutiveScalar*>(itOutput.second.get());
            tangent(0, 0) = mCt * mRho;
            break;
        }
        case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            // nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

bool HeatConduction::CheckHaveParameter(Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
    case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}

double HeatConduction::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
        return this->mK;
    case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
        return this->mCt;
    case Constitutive::eConstitutiveParameter::DENSITY:
        return this->mRho;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void HeatConduction::SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    ConstitutiveBase::CheckParameterDouble(rIdentifier, rValue);
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY:
        this->mK = rValue;
        break;
    case Constitutive::eConstitutiveParameter::DENSITY:
        this->mRho = rValue;
        break;
    case Constitutive::eConstitutiveParameter::HEAT_CAPACITY:
        this->mCt = rValue;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

bool HeatConduction::CheckOutputTypeCompatibility(Constitutive::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::eOutput::HEAT_FLUX:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}

Constitutive::eConstitutiveType HeatConduction::GetType() const
{
    return Constitutive::eConstitutiveType::HEAT_CONDUCTION;
}


void HeatConduction::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Thermal conductivity          : " << this->mK << "\n";
    rLogger << "    Heat capacity                 : " << this->mCt << "\n";
    rLogger << "    Density                       : " << this->mRho << "\n";
}

void HeatConduction::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, mK);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::HEAT_CAPACITY, mCt);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::HEAT_CAPACITY, mRho);
}

template void HeatConduction::Evaluate<1>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void HeatConduction::Evaluate<2>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void HeatConduction::Evaluate<3>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
