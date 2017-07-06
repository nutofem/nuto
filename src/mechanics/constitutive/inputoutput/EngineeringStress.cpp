#include "mechanics/constitutive/inputoutput/EngineeringStress.h"


namespace NuTo
{

template <>
EngineeringStress<3> EngineeringStress<1>::As3D(ePlaneState) const
{
    EngineeringStress<3> stress;
    stress[0] = (*this)[0];
    return stress;
}

template <>
EngineeringStress<3> EngineeringStress<2>::As3D(ePlaneState rPlaneState) const
{
    EngineeringStress<3> stress;
    stress[0] = (*this)[0];
    stress[1] = (*this)[1];
    //  stress[2] = 0.;
    //  stress[3] = 0.;
    //  stress[4] = 0.;
    stress[5] = (*this)[2];
    if (rPlaneState == ePlaneState::PLANE_STRAIN)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 "Not implemented for PLANE_STRAIN and I don't know how to solve it.");
    }
    return stress;
}

template <>
EngineeringStress<3> EngineeringStress<3>::As3D(ePlaneState) const
{
    return *this;
}

template <>
double EngineeringStress<1>::GetVonMisesStress(ePlaneState) const
{
    return (*this)[0];
}

template <>
double EngineeringStress<2>::GetVonMisesStress(ePlaneState rPlaneState) const
{
    if (rPlaneState == ePlaneState::PLANE_STRAIN)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 "Not implemented for PLANE_STRAIN and I don't know how to solve it.");
    }
    const auto& s = data();
    double misesSquared = s[0] * s[0] - s[0] * s[1] + s[1] * s[1] + 3 * s[2] * s[2];
    return std::sqrt(misesSquared);
}


template <>
double EngineeringStress<3>::GetVonMisesStress(ePlaneState) const
{
    const auto& s = data();
    double misesSquared = 0.5 * ((s[0] - s[1]) * (s[0] - s[1]) + (s[1] - s[2]) * (s[1] - s[2]) +
                                 (s[2] - s[0]) * (s[2] - s[0]) + 6 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
    return std::sqrt(misesSquared);
}


} /* namespace NuTo */


template class NuTo::EngineeringStress<1>;
template class NuTo::EngineeringStress<2>;
template class NuTo::EngineeringStress<3>;
