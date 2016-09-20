#include "nuto/mechanics/constitutive/laws/ThermalStrains.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

using namespace NuTo;

template <int TDim>
Error::eError ThermalStrains::Evaluate(const ConstitutiveInputMap &rConstitutiveInput,
        const ConstitutiveOutputMap &rConstitutiveOutput, Constitutive::StaticData::Component*)
{
    
    auto eye = Eigen::MatrixXd::Identity(TDim, TDim);
    const int voigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    double temperature = 0.0;
    std::array<double, 2> strain;
    try
    {
        temperature = (*rConstitutiveInput.at(Constitutive::Input::TEMPERATURE))[0];
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
        case Constitutive::Output::ENGINEERING_STRAIN:
        {
            Eigen::Matrix<double, voigtDim, 1>& engineeringStrain =
                static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second.get())->AsVector();
            for(unsigned int i = 0; i < TDim; ++i)
                engineeringStrain[i] = strain[0];
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case Constitutive::Output::THERMAL_STRAIN:
        {
            Eigen::Matrix<double, TDim, TDim>& engineeringStrain =
                static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second.get())->AsMatrix();
            engineeringStrain = strain[0] * eye;
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case Constitutive::Output::D_STRAIN_D_TEMPERATURE:
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

        case Constitutive::Output::UPDATE_STATIC_DATA:
        case Constitutive::Output::UPDATE_TMP_STATIC_DATA:
            break;
        default:
            continue;
        }
    }
    return Error::SUCCESSFUL;
}

bool ThermalStrains::CheckDofCombinationComputable(Node::eDof rDofRow,
            Node::eDof rDofCol, int) const
{
    if(Node::CombineDofs(rDofRow, rDofCol) == Node::CombineDofs(Node::DISPLACEMENTS, Node::TEMPERATURE))
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

        case Constitutive::Output::ENGINEERING_STRAIN:
        case Constitutive::Output::THERMAL_STRAIN:
        case Constitutive::Output::D_STRAIN_D_TEMPERATURE:
            constitutiveInputMap[Constitutive::Input::TEMPERATURE];
            break;

        case Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case Constitutive::Output::UPDATE_STATIC_DATA:
            break;
        default:
            break;
        }
    }

    return constitutiveInputMap;
}

bool ThermalStrains::CheckElementCompatibility(Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case Element::CONTINUUMELEMENT:
    case Element::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
        return true;
    default:
        return false;
    }
}

void ThermalStrains::SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue)
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

void ThermalStrains::SetParameterFunction(std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    NonlinearExpansionCoeff = ExpansionFunction;
}
