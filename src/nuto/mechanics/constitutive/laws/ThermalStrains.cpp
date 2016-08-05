#include "nuto/mechanics/constitutive/laws/ThermalStrains.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

template <int TDim>
NuTo::Error::eError NuTo::ThermalStrains::Evaluate(NuTo::ElementBase *rElement,
        int rIntegrationPoint, const NuTo::ConstitutiveInputMap &rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    
    auto eye = Eigen::MatrixXd::Identity(TDim, TDim);
    const int voigtDim = NuTo::ConstitutiveIOBase::GetVoigtDim(TDim);

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

    for (auto itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN:
        {
            Eigen::Matrix<double, voigtDim, 1>& engineeringStrain =
                (*static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second)).AsVector();
            for(unsigned int i = 0; i < TDim; ++i)
                engineeringStrain[i] += strain[0];
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case NuTo::Constitutive::Output::THERMAL_STRAIN:
        {
            Eigen::Matrix<double, TDim, TDim>& engineeringStrain =
                (*static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second));
            engineeringStrain = strain[0] * eye;
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case NuTo::Constitutive::Output::D_STRAIN_D_TEMPERATURE:
        {
            Eigen::Matrix<double, voigtDim, 1>& dStrainDTemperature =
                (*static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second)).AsVector();
            for(unsigned int i=0; i<TDim; ++i)
            {
                //! \todo derivative is really positive, yet the strain itself
                //! needs to be negative in the corresponding hessian
                dStrainDTemperature[i] = -strain[1];
            }
            itOutput.second->SetIsCalculated(true);
        }
            break;

        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
            break;
        default:
            continue;
        }
    }
    return NuTo::Error::SUCCESSFUL;
}

bool NuTo::ThermalStrains::CheckDofCombinationComputable(Node::eDof rDofRow,
            Node::eDof rDofCol, int rTimeDerivative) const
{
    if(Node::CombineDofs(rDofRow, rDofCol) == Node::CombineDofs(Node::DISPLACEMENTS, Node::TEMPERATURE))
        return true;

    return false;
}

NuTo::ConstitutiveInputMap NuTo::ThermalStrains::GetConstitutiveInputs(
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
        const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {

        case NuTo::Constitutive::Output::ENGINEERING_STRAIN:
        case NuTo::Constitutive::Output::THERMAL_STRAIN:
        case NuTo::Constitutive::Output::D_STRAIN_D_TEMPERATURE:
            constitutiveInputMap[Constitutive::Input::TEMPERATURE];
            break;

        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
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
    case NuTo::Element::CONTINUUMELEMENT:
    case NuTo::Element::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
        return true;
    default:
        return false;
    }
}

void NuTo::ThermalStrains::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
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

void NuTo::ThermalStrains::SetParameterFunction(std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    NonlinearExpansionCoeff = ExpansionFunction;
}
