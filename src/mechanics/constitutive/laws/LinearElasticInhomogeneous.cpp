#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/LinearElasticInhomogeneous.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include <stdlib.h>
#include <iostream>

void NuTo::LinearElasticInhomogeneous::UpdateParameters(Eigen::VectorXd coordinates)
{
    SetParameterDouble(
                Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,mE(coordinates));
}

void NuTo::LinearElasticInhomogeneous::SetYoungsModulus(std::function<double(Eigen::VectorXd)> E)
{
    mE = E;
}

NuTo::ConstitutiveInputMap NuTo::LinearElasticInhomogeneous::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
            constitutiveInputMap[Constitutive::eInput::COORDINATES];
            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
            constitutiveInputMap[Constitutive::eInput::COORDINATES];
            break;
        }
        // no inputs needed for:
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            constitutiveInputMap[Constitutive::eInput::COORDINATES];
            break;
        }
        default:
            continue;
        }
    }

    return constitutiveInputMap;
}

namespace NuTo // template specialization in same namespace as definition
{

template <int TDim>
void LinearElasticInhomogeneous::Evaluate(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    Eigen::MatrixXd coordinates =
        rConstitutiveInput.at(Constitutive::eInput::COORDINATES)->CopyToEigenMatrix();
    UpdateParameters(coordinates);

    LinearElasticEngineeringStress::Evaluate<TDim>(rConstitutiveInput, rConstitutiveOutput);
}

template void LinearElasticInhomogeneous::Evaluate<1>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void LinearElasticInhomogeneous::Evaluate<2>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void LinearElasticInhomogeneous::Evaluate<3>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);

} // namespace NuTo
