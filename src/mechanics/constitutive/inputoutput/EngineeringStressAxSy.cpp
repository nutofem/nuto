#include "mechanics/constitutive/inputoutput/EngineeringStressAxSy.h"
#include <eigen3/Eigen/Eigenvalues> 

NuTo::EngineeringStressAxSy::EngineeringStressAxSy(std::initializer_list<double> initList)
    : ConstitutiveVector<4>(initList)
{
}

Eigen::Matrix3d ToTensor(const NuTo::EngineeringStressAxSy& stress)
{
    Eigen::Matrix3d s;
    s.setZero();
    s(0, 0) = stress[0];
    s(1, 1) = stress[1];
    s(2, 2) = stress[2];

    s(1, 0) = stress[3];
    s(0, 1) = stress[3];
    return s;
}

//Eigen::Matrix3d ToTensor(const NuTo::EngineeringStress<3>& stress)
//{
//    Eigen::Matrix3d s;
//    s.setZero();
//    s(0, 0) = stress[0];
//    s(1, 1) = stress[1];
//    s(2, 2) = stress[2];
//
//    s(1, 0) = stress[5];
//    s(0, 1) = stress[5];
//
//    s(2, 0) = stress[4];
//    s(0, 2) = stress[4];
//
//    s(2, 1) = stress[3];
//    s(1, 2) = stress[3];
//    return s;
//}

double NuTo::EngineeringStressAxSy::SmoothRankine() const
{
//    Eigen::Matrix3d stressTensor = ToTensor(As3D());
    Eigen::Matrix3d stressTensor = ToTensor(*this);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver;
    Eigen::Vector3d eigenvalues = eigensolver.compute(stressTensor).eigenvalues();
    auto positiveEigenValues = eigenvalues.cwiseMax(Eigen::Vector3d::Zero());
    return positiveEigenValues.norm();
}

namespace NuTo
{

EngineeringStress<3> EngineeringStressAxSy::As3D() const
{
	EngineeringStress<3> stress;

	stress[0] = (*this)[0];
	stress[1] = (*this)[1];
	stress[2] = (*this)[2];
	stress[3] = 0.;
	stress[4] = 0.;
	stress[5] = (*this)[3];

    return stress;
}



double EngineeringStressAxSy::VonMisesStress() const
{
    const auto& s = data();
    double misesSquared = 0.5 * ((s[0] - s[1]) * (s[0] - s[1]) + (s[1] - s[2]) * (s[1] - s[2]) +
                                 (s[2] - s[0]) * (s[2] - s[0]) + 6 * s[3] * s[3]);
    return std::sqrt(misesSquared);
}

} /* namespace NuTo */
