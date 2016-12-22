#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "base/ErrorEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"


NuTo::Constitutive::eConstitutiveType NuTo::AdditiveOutput::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT;
}

template NuTo::eError NuTo::AdditiveOutput::Evaluate<1>(const NuTo::ConstitutiveInputMap &rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput);
template NuTo::eError NuTo::AdditiveOutput::Evaluate<2>(const NuTo::ConstitutiveInputMap &rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput);
template NuTo::eError NuTo::AdditiveOutput::Evaluate<3>(const NuTo::ConstitutiveInputMap &rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput);
