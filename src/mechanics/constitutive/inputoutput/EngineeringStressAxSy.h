#pragma once

#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"

namespace NuTo
{

//! @brief ... Engineering stress
/*!
 *  2D case, axisymmetric:
 *  \f[
 *     \boldsymbol{\sigma} = \begin{bmatrix}
 *        \sigma_{r}\\
 *        \sigma_{z}\\
 *        \sigma_{th}\\
 *        \tau_{rz}
 *    \end{bmatrix} = \begin{bmatrix}
 *        \sigma_{rr}\\
 *        \sigma_{zz}\\
 *        \sigma_{thth}\\
 *        \sigma_{rz}
 *    \end{bmatrix}.
 *  \f]
 */
class EngineeringStressAxSy : public ConstitutiveVector<4>
{
public:
    EngineeringStressAxSy() = default;
    EngineeringStressAxSy(std::initializer_list<double> initList);
    //EngineeringStressAxSy(const EngineeringStressAxSy&) = default;
    //EngineeringStressAxSy(EngineeringStressAxSy&&) = default;

    EngineeringStress<3> As3D() const;

    double VonMisesStress() const;

    double SmoothRankine() const;
};

} /* namespace NuTo */
