#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"


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
