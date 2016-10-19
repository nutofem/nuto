#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/laws/ThermalStrains.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

using namespace NuTo;

template <int TDim>
NuTo::eError ThermalStrains::Evaluate(
        const ConstitutiveInputMap &rConstitutiveInput,
        const ConstitutiveOutputMap &rConstitutiveOutput,
        Constitutive::StaticData::Component*)
{
    
    auto eye = Eigen::MatrixXd::Identity(TDim, TDim);
    const int voigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    double temperature = 0.0;
    std::array<double, 2> strain;
    try
    {
        temperature = (*rConstitutiveInput.at(Constitutive::eInput::TEMPERATURE))[0];
        strain = NonlinearExpansionCoeff(temperature);
    }
    catch(const std::bad_function_call& e)
    {
        strain[0] = mExpansionCoefficient * temperature;
    }
    catch(const std::out_of_range& e) {}

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
            Eigen::Matrix<double, TDim, TDim>& engineeringStrain =
                static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second.get())->AsMatrix();
            engineeringStrain = strain[0] * eye;
            itOutput.second->SetIsCalculated(true);
            break;
        }
        case NuTo::Constitutive::eOutput::D_STRAIN_D_TEMPERATURE:
        {
            Eigen::Matrix<double, voigtDim, 1>& dStrainDTemperature =
                static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second.get())->AsVector();
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

bool ThermalStrains::CheckDofCombinationComputable(Node::eDof rDofRow,
            Node::eDof rDofCol, int) const
{
    if(Node::CombineDofs(rDofRow, rDofCol) == Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::TEMPERATURE))
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
    NonlinearExpansionCoeff = ExpansionFunction;
}
