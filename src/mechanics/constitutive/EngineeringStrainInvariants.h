#pragma once

#include <eigen3/Eigen/Core>
#include "mechanics/constitutive/EngineeringStrain.h"

namespace NuTo
{
namespace EngineeringStrainInvariants
{
using Tensor = Eigen::Matrix3d;

Tensor ToTensor(EngineeringStrain<3>& v)
{
    Tensor t = Tensor::Zero();
    t(0, 0) = v[0];
    t(1, 1) = v[1];
    t(2, 2) = v[2];

    t(1, 0) = 0.5 * v[5];
    t(0, 1) = 0.5 * v[5];

    t(2, 0) = 0.5 * v[4];
    t(0, 2) = 0.5 * v[4];

    t(2, 1) = 0.5 * v[3];
    t(1, 2) = 0.5 * v[3];
    return t;
}

//! @brief returns I1 - the first strain invariant of the characteristic equation
//! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
double I1(const EngineeringStrain<3>& v)
{
    return v.segment<3>(0).sum();
}

//! @brief returns I2 - the first strain invariant of the characteristic equation
//! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
double I2(const EngineeringStrain<3>& v)
{
    // e_xx = v[0]     || e_yy = v[1]     || e_zz = v[2]
    // e_yz = v[3] / 2 || e_zx = v[4] / 2 || e_xy = v[5] / 2
    return v[0] * v[1] + v[1] * v[2] + v[2] * v[0] - 0.25 * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);
}

//! @brief returns I3 - the first strain invariant of the characteristic equation
//! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
double I3(const EngineeringStrain<3>& v)
{
    // e_xx = v[0]     || e_yy = v[1]     || e_zz = v[2]
    // e_yz = v[3] / 2 || e_zx = v[4] / 2 || e_xy = v[5] / 2
    return v[0] * v[1] * v[2] + 0.25 * v[3] * v[4] * v[5] -
           .25 * (v[5] * v[5] * v[2] + v[3] * v[3] * v[0] + v[4] * v[4] * v[1]);
}

//! @brief returns J2 - the second deviatoric strain invariant of the characteristic equation
//! \f[ \lambda^3 - J_1 \lambda^2 - J_2 \lambda - J_3  \f] Note the minus sign in front of J2. This is not
//! consistent with the I2 invariant but apparently common practice.
double J2(const EngineeringStrain<3>& v)
{
    const double eps_xx = v[0];
    const double eps_yy = v[1];
    const double eps_zz = v[2];
    const double eps_yz = .5 * v[3];
    const double eps_zx = .5 * v[4];
    const double eps_xy = .5 * v[5];

    return 1. / 6. * ((eps_xx - eps_yy) * (eps_xx - eps_yy) + (eps_yy - eps_zz) * (eps_yy - eps_zz) +
                      (eps_zz - eps_xx) * (eps_zz - eps_xx)) +
           eps_xy * eps_xy + eps_yz * eps_yz + eps_zx * eps_zx;
}


//! @brief returns the deviatoric part
EngineeringStrain<3> Deviatoric(EngineeringStrain<3> v)
{
    double I1_3 = I1(v) / 3.;
    v[0] -= I1_3;
    v[1] -= I1_3;
    v[2] -= I1_3;
    return v;
}

} /* TensorCalc */
} /* NuTo */
