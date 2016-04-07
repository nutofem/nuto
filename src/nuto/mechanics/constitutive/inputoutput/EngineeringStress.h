#pragma once

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"


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
class EngineeringStress: public ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>
{
public:

    double GetVonMisesStress() const;

};

} /* namespace NuTo */

