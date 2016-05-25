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

    for (auto itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN:
        {
            double temperature = (*rConstitutiveInput.at(Constitutive::Input::TEMPERATURE))[0];
            Eigen::Matrix<double, voigtDim, 1>& engineeringStrain =
                (*static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second)).AsVector();
            for(unsigned int i = 0; i < TDim; ++i)
                engineeringStrain[i] += mExpansionCoefficient * temperature;
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case NuTo::Constitutive::Output::THERMAL_STRAIN_VISUALIZE:
        {
            double temperature = (*rConstitutiveInput.at(Constitutive::Input::TEMPERATURE))[0];
            Eigen::Matrix<double, TDim, TDim>& engineeringStrain =
                (*static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second));
            engineeringStrain = mExpansionCoefficient * eye * temperature;
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case NuTo::Constitutive::Output::D_STRAIN_D_TEMPERATURE:
        {
            Eigen::Matrix<double, voigtDim, 1>& dStrainDTemperature =
                (*static_cast<ConstitutiveVector<voigtDim>*>(itOutput.second)).AsVector();
            for(unsigned int i=0; i<TDim; ++i)
            {
                //TODO: This can't be right; need to investigate
                //dStrainDTemperature[i]+= mExpansionCoefficient;
                dStrainDTemperature[i] = 0.0;
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

