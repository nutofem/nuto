#include "LinearDampingEngineeringStress.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawWithoutData.h"


NuTo::LinearDampingEngineeringStress::LinearDampingEngineeringStress()
    : ConstitutiveBase()
    , mDampingCoefficient(-1.0)
{
}

std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> NuTo::LinearDampingEngineeringStress::CreateIPLaw()
{
    return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<LinearDampingEngineeringStress>>(*this);
}


NuTo::ConstitutiveInputMap
NuTo::LinearDampingEngineeringStress::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                            const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN_DT1];
            break;

        // no inputs needed for:
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            continue;
        }
    }

    return constitutiveInputMap;
}


bool NuTo::LinearDampingEngineeringStress::CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
                                                                         int rTimeDerivative) const
{
    assert(rTimeDerivative == 0 or rTimeDerivative == 1 or rTimeDerivative == 2);
    switch (rTimeDerivative)
    {
    case 1:
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


void NuTo::LinearDampingEngineeringStress::CheckParameters() const
{
}


NuTo::Constitutive::eConstitutiveType NuTo::LinearDampingEngineeringStress::GetType() const
{
    return Constitutive::eConstitutiveType::LINEAR_DAMPING_ENGINEERING_STRESS;
}


bool NuTo::LinearDampingEngineeringStress::HaveTmpStaticData() const
{
    return false;
}


double NuTo::LinearDampingEngineeringStress::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT:
        return mDampingCoefficient;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the parameter " +
                                                              Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}


void NuTo::LinearDampingEngineeringStress::SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier,
                                                              double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT:
        mDampingCoefficient = rValue;
        break;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the parameter " +
                                                              Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
    SetParametersValid();
}


template <int TDim>
void NuTo::LinearDampingEngineeringStress::Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                                                    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            assert(rConstitutiveInput.find(Constitutive::eInput::ENGINEERING_STRAIN_DT1) != rConstitutiveInput.end());
            const auto& engineeringStrainDT1 =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN_DT1)->AsEngineeringStrain<TDim>();

            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<TDim>(itOutput.first, __PRETTY_FUNCTION__);
            for (unsigned int i = 0; i < TDim; ++i)
            {
                engineeringStress[i] = mDampingCoefficient * engineeringStrainDT1[i];
            }
            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1:
        {
            if (TDim == 1)
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<1, 1>(itOutput.first, __PRETTY_FUNCTION__);

                tangent(0, 0) = mDampingCoefficient;
            }
            else
                throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Not implemented for dimension 2 and 3");

            break;
        }

        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            continue;

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

template void NuTo::LinearDampingEngineeringStress::Evaluate<1>(const ConstitutiveInputMap& rConstitutiveInput,
                                                                const ConstitutiveOutputMap& rConstitutiveOutput);
template void NuTo::LinearDampingEngineeringStress::Evaluate<2>(const ConstitutiveInputMap& rConstitutiveInput,
                                                                const ConstitutiveOutputMap& rConstitutiveOutput);
template void NuTo::LinearDampingEngineeringStress::Evaluate<3>(const ConstitutiveInputMap& rConstitutiveInput,
                                                                const ConstitutiveOutputMap& rConstitutiveOutput);
