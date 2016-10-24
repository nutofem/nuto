#include "nuto/mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

std::tuple<double, double, double> NuTo::EngineeringStressHelper::CalculateCoefficients2DPlaneStress(double rE, double rNu)
{
    double factor = rE / (1.0 - (rNu * rNu));
    return std::make_tuple(
            factor,                         // C11
            factor * rNu,                   // C12
            factor * 0.5 * (1.0 - rNu));    // C33
}


std::tuple<double, double, double> NuTo::EngineeringStressHelper::CalculateCoefficients3D(double rE, double rNu)
{
    double factor = rE / ((1.0 + rNu) * (1.0 - 2.0 * rNu));
    return std::make_tuple(
            factor * (1.0 - rNu),           // C11
            factor * rNu,                   // C12
            rE / (2. * (1.0 + rNu)));       // C33
}


namespace NuTo
{
template <>
EngineeringStress<1> EngineeringStressHelper::GetStress<1>(const EngineeringStrain<1>& rElasticStrain,
                                                           double rE,
                                                           double rNu, ePlaneState rPlaneState)
{
    EngineeringStress<1> stress;
    stress[0] = rE * rElasticStrain[0];
    return stress;
}

template <>
EngineeringStress<2> EngineeringStressHelper::GetStress<2>(const EngineeringStrain<2>& rElasticStrain,
                                                           double rE,
                                                           double rNu, ePlaneState rPlaneState)
{
    EngineeringStress<2> stress;
    if (rPlaneState == ePlaneState::PLANE_STRESS)
    {

    }
    else
    {

    }
    return stress;
}

template <>
EngineeringStress<3> EngineeringStressHelper::GetStress<3>(const EngineeringStrain<3>& rElasticStrain,
                                                           double rE,
                                                           double rNu, ePlaneState rPlaneState)
{
    EngineeringStress<3> stress;
    return stress;
}


} // namespace NuTo