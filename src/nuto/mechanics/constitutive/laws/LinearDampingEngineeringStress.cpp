#include "LinearDampingEngineeringStress.h"

#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

NuTo::LinearDampingEngineeringStress::LinearDampingEngineeringStress() : ConstitutiveBase() {}


NuTo::ConstitutiveInputMap NuTo::LinearDampingEngineeringStress::GetConstitutiveInputs(
        const ConstitutiveOutputMap &rConstitutiveOutput, const InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
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


bool NuTo::LinearDampingEngineeringStress::CheckDofCombinationComputable(Node::eDof rDofRow,
                                                                         Node::eDof rDofCol,
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


bool NuTo::LinearDampingEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
    case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
    case NuTo::Element::eElementType::ELEMENT1DINXD:
        return true;
    default:
        return false;
    }
}


void NuTo::LinearDampingEngineeringStress::CheckParameters() const {}


NuTo::eError NuTo::LinearDampingEngineeringStress::Evaluate1D(const ConstitutiveInputMap &rConstitutiveInput,
                                                              const ConstitutiveOutputMap &rConstitutiveOutput,
                                                              Constitutive::StaticData::Component* staticData)
{
    return EvaluateLinearDampingEngineeringStress<1>(rConstitutiveInput, rConstitutiveOutput, staticData);
}


NuTo::eError NuTo::LinearDampingEngineeringStress::Evaluate2D(const ConstitutiveInputMap &rConstitutiveInput,
                                                              const ConstitutiveOutputMap &rConstitutiveOutput,
                                                              Constitutive::StaticData::Component* staticData)
{
    return EvaluateLinearDampingEngineeringStress<2>(rConstitutiveInput, rConstitutiveOutput, staticData);
}


NuTo::eError NuTo::LinearDampingEngineeringStress::Evaluate3D(const ConstitutiveInputMap &rConstitutiveInput,
                                                              const ConstitutiveOutputMap &rConstitutiveOutput,
                                                              Constitutive::StaticData::Component* staticData)
{
    return EvaluateLinearDampingEngineeringStress<3>(rConstitutiveInput, rConstitutiveOutput, staticData);
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
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT:
        return mDampingCoefficient;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,
                "Constitutive law does not have the parameter "
                + Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}


void NuTo::LinearDampingEngineeringStress::SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier,
        double rValue)
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT:
        mDampingCoefficient = rValue;
        break;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,
                "Constitutive law does not have the parameter "
                + Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
    SetParametersValid();
}


template <int TDim>
NuTo::eError NuTo::LinearDampingEngineeringStress::EvaluateLinearDampingEngineeringStress(
        const ConstitutiveInputMap &rConstitutiveInput,
        const ConstitutiveOutputMap &rConstitutiveOutput,
        Constitutive::StaticData::Component*)
{
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            const auto& engineeringStrainDT1 =
                rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN_DT1)->AsEngineeringStrain<TDim>();

            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<TDim>(itOutput.first, __PRETTY_FUNCTION__);
            for(unsigned int i=0; i<TDim; ++i)
            {
                engineeringStress[i] = mDampingCoefficient * engineeringStrainDT1[i];
            }
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

    return eError::SUCCESSFUL;
}
