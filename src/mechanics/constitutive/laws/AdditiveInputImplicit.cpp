#include "AdditiveInputImplicit.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

#include "mechanics/nodes/NodeEnum.h"

#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseDirectSolverMKLPardiso.h"
#include "math/SparseDirectSolverPardiso.h"
#include "math/SparseMatrixCSRGeneral.h"


NuTo::Constitutive::eConstitutiveType NuTo::AdditiveInputImplicit::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_IMPLICIT;
}

NuTo::ConstitutiveInputMap
NuTo::AdditiveInputImplicit::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap& rConstitutiveOutput,
                                                   const NuTo::InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap =
            AdditiveBase::GetConstitutiveInputs(rConstitutiveOutput, rInterpolationType);

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::ENGINEERING_STRESS:
        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            break;
        default:
            break;
        }
    }
    return constitutiveInputMap;
}


template void NuTo::AdditiveInputImplicit::Evaluate<1>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                       const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);

template void NuTo::AdditiveInputImplicit::Evaluate<2>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                       const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);

template void NuTo::AdditiveInputImplicit::Evaluate<3>(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                       const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
