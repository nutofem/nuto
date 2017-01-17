#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

NuTo::Constitutive::eConstitutiveType NuTo::AdditiveOutput::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT;
}
