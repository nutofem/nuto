#include "AdditiveInputImplicit.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"
#include "nuto/math/SparseDirectSolverPardiso.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"




NuTo::Constitutive::eConstitutiveType NuTo::AdditiveInputImplicit::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_IMPLICIT;
}

NuTo::ConstitutiveInputMap NuTo::AdditiveInputImplicit::GetConstitutiveInputs(
        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput, const NuTo::InterpolationType &rInterpolationType) const
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


template NuTo::eError NuTo::AdditiveInputImplicit::Evaluate<1>(
    const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput);

template NuTo::eError NuTo::AdditiveInputImplicit::Evaluate<2>(
    const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput);

template NuTo::eError NuTo::AdditiveInputImplicit::Evaluate<3>(
    const NuTo::ConstitutiveInputMap &rConstitutiveInput, const NuTo::ConstitutiveOutputMap &rConstitutiveOutput);
