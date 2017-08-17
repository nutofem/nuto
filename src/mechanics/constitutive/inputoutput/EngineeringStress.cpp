#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include <eigen3/Eigen/Eigenvalues> 

template <int TDim>
NuTo::EngineeringStress<TDim>::EngineeringStress(std::initializer_list<double> initList)
    : ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>(initList)
{
}

Eigen::Matrix3d ToTensor(const NuTo::EngineeringStress<3>& stress)
{
    Eigen::Matrix3d s;
    s.setZero();
    s(0, 0) = stress[0];
    s(1, 1) = stress[1];
    s(2, 2) = stress[2];

    s(1, 0) = stress[5];
    s(0, 1) = stress[5];

    s(2, 0) = stress[4];
    s(0, 2) = stress[4];

    s(2, 1) = stress[3];
    s(1, 2) = stress[3];
    return s;
}

template <int TDim>
double NuTo::EngineeringStress<TDim>::SmoothRankine(ePlaneState planeState) const
{
    Eigen::Matrix3d stressTensor = ToTensor(As3D());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver;
    Eigen::Vector3d eigenvalues = eigensolver.compute(stressTensor).eigenvalues();
    auto positiveEigenValues = eigenvalues.cwiseMax(Eigen::Vector3d::Zero());
    return positiveEigenValues.norm();
}

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
        throw Exception(__PRETTY_FUNCTION__, "Not implemented for PLANE_STRAIN and I don't know how to solve it.");
    }
    return stress;
}

template <>
EngineeringStress<3> EngineeringStress<3>::As3D(ePlaneState) const
{
    return *this;
}

template <>
double EngineeringStress<1>::VonMisesStress(ePlaneState) const
{
    return (*this)[0];
}

template <>
double EngineeringStress<2>::VonMisesStress(ePlaneState rPlaneState) const
{
    if (rPlaneState == ePlaneState::PLANE_STRAIN)
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented for PLANE_STRAIN and I don't know how to solve it.");
    }
    const auto& s = data();
    double misesSquared = s[0] * s[0] - s[0] * s[1] + s[1] * s[1] + 3 * s[2] * s[2];
    return std::sqrt(misesSquared);
}


template <>
double EngineeringStress<3>::VonMisesStress(ePlaneState) const
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
