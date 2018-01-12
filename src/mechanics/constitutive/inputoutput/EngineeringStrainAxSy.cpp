#include "mechanics/constitutive/inputoutput/EngineeringStrainAxSy.h"

NuTo::EngineeringStrainAxSy::EngineeringStrainAxSy(std::initializer_list<double> initList)
    : ConstitutiveVector<4>(initList)
{
}


//std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::EngineeringStrainAxSy::clone()
//{
//    return std::make_unique<EngineeringStrainAxSy>(*this);
//}

double NuTo::EngineeringStrainAxSy::InvariantI1() const
{
	throw Exception(__PRETTY_FUNCTION__, "Not implemented for axisymmetric");
//    return (*this).segment<3>(0).sum(); // sum of the first 3 entries
}

double NuTo::EngineeringStrainAxSy::InvariantI2() const
{
	throw Exception(__PRETTY_FUNCTION__, "Not implemented for axisymmetric");

//    const Eigen::Matrix<double, 4, 1>& e = *this;
//    // e_rr   = v[0]
//    // e_zz   = v[1]
//    // e_thth = v[2]
//    // e_rz   = v[3] / 2
//    return e[0] * e[1] + e[1] * e[2] + e[2] * e[0] - 0.25 * e[3] * e[3];
}

double NuTo::EngineeringStrainAxSy::InvariantI3() const
{
	throw Exception(__PRETTY_FUNCTION__, "Not implemented for axisymmetric");

//    const Eigen::Matrix<double, 4, 1>& v = *this;
    // e_xx = v[0]
    // e_yy = v[1]
    // e_zz = v[2]
    // e_yz = v[3] / 2
    // e_zx = v[4] / 2
    // e_xy = v[5] / 2
//    // e_rr   = v[0]
//    // e_zz   = v[1]
//    // e_thth = v[2]
//    // e_rz   = v[3] / 2
//    return v[0] * v[1] * v[2] - .25 * v[3] * v[3] * v[0];
}

double NuTo::EngineeringStrainAxSy::InvariantJ2() const
{
	throw Exception(__PRETTY_FUNCTION__, "Not implemented for axisymmetric");

//    const Eigen::Matrix<double, 4, 1>& v = *this;
//    const double eps_xx = v[0];
//    const double eps_yy = v[1];
//    const double eps_zz = v[2];
//    const double eps_yz = .5 * v[3];
//    const double eps_zx = .5 * 0.;
//    const double eps_xy = .5 * 0.;
//
//    return 1. / 6. * ((eps_xx - eps_yy) * (eps_xx - eps_yy) + (eps_yy - eps_zz) * (eps_yy - eps_zz) +
//                      (eps_zz - eps_xx) * (eps_zz - eps_xx)) +
//           eps_xy * eps_xy + eps_yz * eps_yz + eps_zx * eps_zx;
}


NuTo::EngineeringStrainAxSy NuTo::EngineeringStrainAxSy::Deviatoric() const
{
	throw Exception(__PRETTY_FUNCTION__, "Not implemented for axisymmetric");

//    EngineeringStrainAxSy v = *this;
//    double I1_3 = v.InvariantI1() / 3.;
//    v[0] -= I1_3;
//    v[1] -= I1_3;
//    v[2] -= I1_3;
//    return v;
}

namespace NuTo
{


NuTo::EngineeringStrain<3> NuTo::EngineeringStrainAxSy::As3D() const
{
    EngineeringStrain<3> strain3D;
    const EngineeringStrainAxSy& v = *this;

    strain3D[0] = v[0];
    strain3D[1] = v[1];
    strain3D[2] = v[2];
    strain3D[3] = 0.;
    strain3D[4] = 0.;
    strain3D[5] = v[3];

    return strain3D;
}
} // namespace NuTo
