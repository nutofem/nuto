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
    double C11 = 0., C12 = 0., C33 = 0.;
    if (rPlaneState == ePlaneState::PLANE_STRESS)
        std::tie(C11, C12, C33) = CalculateCoefficients3D(rE, rNu);
    else
        std::tie(C11, C12, C33) = CalculateCoefficients2DPlaneStress(rE, rNu);

    EngineeringStress<2> stress;
    stress[0] = C11 * rElasticStrain[0] + C12 * rElasticStrain[1];
    stress[1] = C11 * rElasticStrain[1] + C12 * rElasticStrain[0];
    stress[2] = C33 * rElasticStrain[2];
    return stress;
}

template <>
EngineeringStress<3> EngineeringStressHelper::GetStress<3>(const EngineeringStrain<3>& rElasticStrain,
                                                           double rE,
                                                           double rNu, ePlaneState rPlaneState)
{
    double C11 = 0.0, C12 = 0.0, C44 = 0.0;
    std::tie(C11, C12, C44) = EngineeringStressHelper::CalculateCoefficients3D(rE, rNu);

    EngineeringStress<3> stress;
    stress[0] = C11 * rElasticStrain[0] + C12 * (rElasticStrain[1] + rElasticStrain[2]);
    stress[1] = C11 * rElasticStrain[1] + C12 * (rElasticStrain[0] + rElasticStrain[2]);
    stress[2] = C11 * rElasticStrain[2] + C12 * (rElasticStrain[0] + rElasticStrain[1]);
    stress[3] = C44 * rElasticStrain[3];
    stress[4] = C44 * rElasticStrain[4];
    stress[5] = C44 * rElasticStrain[5];
    return stress;
}


} // namespace NuTo