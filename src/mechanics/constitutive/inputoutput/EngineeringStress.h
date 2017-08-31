#pragma once

#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneStateEnum.h"


namespace NuTo
{

//! @brief ... Engineering stress
/*!
 *  3D case:
 *  \f[
 *     \boldsymbol{\sigma} = \begin{bmatrix}
 *        \sigma_{x}\\
 *        \sigma_{y}\\
 *        \sigma_{z}\\
 *        \tau_{yz}\\
 *        \tau_{zx}\\
 *        \tau_{xy}
 *    \end{bmatrix} = \begin{bmatrix}
 *        \sigma_{xx}\\
 *        \sigma_{yy}\\
 *        \sigma_{zz}\\
 *        \sigma_{yz}\\
 *        \sigma_{zx}\\
 *        \sigma_{xy}
 *    \end{bmatrix}.
 *  \f]
 */
template <int TDim>
class EngineeringStress : public ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>
{
public:
    EngineeringStress() = default;
    EngineeringStress(std::initializer_list<double> initList);
    //EngineeringStress(const EngineeringStress&) = default;
    //EngineeringStress(EngineeringStress&&) = default;

    EngineeringStress<3> As3D(ePlaneState rPlaneState = ePlaneState::PLANE_STRESS) const;

    double VonMisesStress(ePlaneState rPlaneState = ePlaneState::PLANE_STRESS) const;

    double SmoothRankine(ePlaneState rPlaneState = ePlaneState::PLANE_STRESS) const;
};

} /* namespace NuTo */
