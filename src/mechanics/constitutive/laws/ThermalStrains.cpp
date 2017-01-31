#include "base/ErrorEnum.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawWithoutData.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> NuTo::ThermalStrains::CreateIPLaw()
{
    return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<ThermalStrains>>(*this);
}


template <int TDim>
NuTo::eError ThermalStrains::Evaluate(
        const ConstitutiveInputMap &rConstitutiveInput,
        const ConstitutiveOutputMap &rConstitutiveOutput)
{
    const int voigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    if (rConstitutiveOutput.size() == 1 and rConstitutiveOutput.count(Constitutive::eOutput::UPDATE_STATIC_DATA) == 1)
        return eError::SUCCESSFUL;

    double temperature = (*rConstitutiveInput.at(Constitutive::eInput::TEMPERATURE))[0];

    std::array<double, 2> strain = {0.0, 0.0};
    if (mNonlinearExpansionFunction)
    {
        if (temperature < 0.0) temperature = 0.0;
        strain = mNonlinearExpansionFunction(temperature);
    }
    else
    {
        strain[0] = mExpansionCoefficient * temperature;
        strain[1] = mExpansionCoefficient;
    }

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN:
        {
            Eigen::Matrix<double, voigtDim, 1>& engineeringStrain =
                static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second.get())->AsVector();
            for(unsigned int i = 0; i < TDim; ++i)
                engineeringStrain[i] = strain[0];
            itOutput.second->SetIsCalculated(true);
            break;
        }
        case NuTo::Constitutive::eOutput::THERMAL_STRAIN:
        {
            ConstitutiveVector<voigtDim>& engineeringStrain =
                *static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second.get());
            for(unsigned int i = 0; i < TDim; ++i)
                engineeringStrain[i] = strain[0];
            itOutput.second->SetIsCalculated(true);
            break;
        }
        case NuTo::Constitutive::eOutput::D_STRAIN_D_TEMPERATURE:
        {
            Eigen::Matrix<double, voigtDim, 1>& dStrainDTemperature =
                static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second.get())->AsVector();
            dStrainDTemperature.setZero();
            for(unsigned int i=0; i<TDim; ++i)
            {
                //! \todo derivative is really positive, yet the strain itself
                //! needs to be negative in the corresponding hessian
                dStrainDTemperature[i] = -strain[1];
            }
            itOutput.second->SetIsCalculated(true);
        }
            break;

        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
            break;
        default:
            continue;
        }
    }
    return NuTo::eError::SUCCESSFUL;
}

bool ThermalStrains::CheckDofCombinationComputable(Node::eDof dofRow, Node::eDof dofCol, int timeDerivative) const
{
    if (timeDerivative > 1)
        return false;

    if (Node::CombineDofs(dofRow, dofCol) == Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE))
        return true;

    return false;
}

ConstitutiveInputMap ThermalStrains::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
        const InterpolationType&) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {

        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN:
        case NuTo::Constitutive::eOutput::THERMAL_STRAIN:
        case NuTo::Constitutive::eOutput::D_STRAIN_D_TEMPERATURE:
            constitutiveInputMap[Constitutive::eInput::TEMPERATURE];
            break;

        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            break;
        }
    }

    return constitutiveInputMap;
}

bool NuTo::ThermalStrains::CheckElementCompatibility(Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
    case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
        return true;
    default:
        return false;
    }
}

void NuTo::ThermalStrains::SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
            mExpansionCoefficient = rValue;
            return;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__,
                    "Constitutive law does not have the parameter "
                    + Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

Constitutive::eConstitutiveType NuTo::ThermalStrains::GetType() const
{
    return Constitutive::eConstitutiveType::THERMAL_STRAINS;
}

void NuTo::ThermalStrains::SetParameterFunction(std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    mNonlinearExpansionFunction = ExpansionFunction;
}

template NuTo::eError ThermalStrains::Evaluate<1>(const ConstitutiveInputMap &rConstitutiveInput,
                                                  const ConstitutiveOutputMap &rConstitutiveOutput);
template NuTo::eError ThermalStrains::Evaluate<2>(const ConstitutiveInputMap &rConstitutiveInput,
                                                  const ConstitutiveOutputMap &rConstitutiveOutput);
template NuTo::eError ThermalStrains::Evaluate<3>(const ConstitutiveInputMap &rConstitutiveInput,
                                                  const ConstitutiveOutputMap &rConstitutiveOutput);