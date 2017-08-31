#pragma once

#include <eigen3/Eigen/Core>
#include "mechanics/constitutive/Voigt.h"

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
class EngineeringStressPDE : public Eigen::Matrix<double, Voigt::Dim(TDim), 1>
{
    using Parent = Eigen::Matrix<double, Voigt::Dim(TDim), 1>;
public:
    EngineeringStressPDE() = default;

    // This constructor allows you to construct EngineeringStressPDE from Eigen expressions
    template <typename OtherDerived>
    EngineeringStressPDE(const Eigen::MatrixBase<OtherDerived>& other)
        : Parent(other) 
    {
    }
    
    // This constructor allows you to construct EngineeringStressPDE from Eigen expressions
    template <typename OtherDerived>
    EngineeringStressPDE(Eigen::MatrixBase<OtherDerived>&& other)
        : Parent(other) 
    {
    }
    
    // This method allows you to assign Eigen expressions to EngineeringStressPDE
    template <typename OtherDerived>
    EngineeringStressPDE& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Parent::operator=(other);
        return *this;
    }
    
    // This method allows you to assign Eigen expressions to EngineeringStressPDE
    template <typename OtherDerived>
    EngineeringStressPDE& operator=(Eigen::MatrixBase<OtherDerived>&& other)
    {
        this->Parent::operator=(other);
        return *this;
    }
};
} /* namespace NuTo */
