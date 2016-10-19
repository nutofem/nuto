#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/laws/HeatConduction.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

NuTo::HeatConduction::HeatConduction() : ConstitutiveBase()
{
    mK = 0.0;
    mCt = 0.0;
    mRho = 0.0;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
template<class Archive>
void NuTo::HeatConduction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize HeatConduction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
    & BOOST_SERIALIZATION_NVP(mK)
    & BOOST_SERIALIZATION_NVP(mCt)
    & BOOST_SERIALIZATION_NVP(mRho);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize HeatConduction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::HeatConduction)
#endif // ENABLE_SERIALIZATION

NuTo::ConstitutiveInputMap NuTo::HeatConduction::GetConstitutiveInputs(
    const ConstitutiveOutputMap& rConstitutiveOutput,
    const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::HEAT_FLUX:
        {
            constitutiveInputMap[Constitutive::eInput::TEMPERATURE_GRADIENT];
            break;
        }
        case NuTo::Constitutive::eOutput::HEAT_CHANGE:
        {
            constitutiveInputMap[Constitutive::eInput::TEMPERATURE_CHANGE];
            break;
        }
        case NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        case NuTo::Constitutive::eOutput::D_HEAT_D_TEMPERATURE:
        {
            // for nonlinear:
            //constitutiveInputMap[Constitutive::eInput::TEMPERATURE];
            break;
        }
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            continue;
        }
    }

    return constitutiveInputMap;
}

bool NuTo::HeatConduction::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    assert(rTimeDerivative>-1);
    if (rTimeDerivative<=2 &&
        rDofRow == Node::eDof::TEMPERATURE &&
        rDofCol == Node::eDof::TEMPERATURE)
    {
        return true;
    }
    return false;
}

template<int TDim>
NuTo::eError NuTo::HeatConduction::Evaluate(
		const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput,
		Constitutive::StaticData::Component*)
{
    auto eye = Eigen::MatrixXd::Identity(TDim, TDim);

    InputData<TDim> inputData;
    for (auto& itInput : rConstitutiveInput)
    {
        switch(itInput.first)
        {
        case NuTo::Constitutive::eInput::TEMPERATURE_GRADIENT:
            inputData.mTemperatureGradient = static_cast<ConstitutiveVector<TDim>*>(itInput.second.get())->AsVector();
            break;
        case NuTo::Constitutive::eInput::TEMPERATURE_CHANGE:
            inputData.mTemperatureChange = (*itInput.second)[0];
            break;
        case NuTo::Constitutive::eInput::CALCULATE_STATIC_DATA:
        case NuTo::Constitutive::eInput::TIME_STEP:
            break;
        default:
            continue;
        }
    }

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::HEAT_FLUX:
        {
            Eigen::Matrix<double, TDim, 1>& heatFlux = *static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get());
            heatFlux = -mK * eye * inputData.mTemperatureGradient;
            break;
        }
        case NuTo::Constitutive::eOutput::HEAT_CHANGE:
        {
            Eigen::Matrix<double, 1, 1>& heatChange = *static_cast<ConstitutiveScalar*>(itOutput.second.get());
            heatChange(0, 0) = mCt * mRho * inputData.mTemperatureChange;
            break;
        }
        case NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        {
            Eigen::Matrix<double, TDim, TDim>& conductivity = *static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second.get());
            conductivity = mK * eye;
            break;
        }
        case NuTo::Constitutive::eOutput::D_HEAT_D_TEMPERATURE:
        {
            Eigen::Matrix<double, 1, 1>& tangent = *static_cast<ConstitutiveScalar*>(itOutput.second.get());
            tangent(0,0) = mCt * mRho;
            break;
        }
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            //nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
    return eError::SUCCESSFUL;
}

bool NuTo::HeatConduction::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
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

double NuTo::HeatConduction::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
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

void NuTo::HeatConduction::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    ConstitutiveBase::CheckParameterDouble(rIdentifier, rValue);
    switch(rIdentifier)
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

bool NuTo::HeatConduction::CheckOutputTypeCompatibility(NuTo::Constitutive::eOutput rOutputEnum) const
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

NuTo::Constitutive::eConstitutiveType NuTo::HeatConduction::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION;
}

bool NuTo::HeatConduction::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
        return true;
    default:
        return false;
    }
}

void NuTo::HeatConduction::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Thermal conductivity          : " << this->mK << "\n";
    rLogger << "    Heat capacity                 : " << this->mCt << "\n";
    rLogger << "    Density                       : " << this->mRho << "\n";
}

void NuTo::HeatConduction::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, mK);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::HEAT_CAPACITY, mCt);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::HEAT_CAPACITY, mRho);
}
