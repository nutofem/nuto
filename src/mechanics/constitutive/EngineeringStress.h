#pragma once

#include <eigen3/Eigen/Core>
#include "mechanics/constitutive/Voigt.h"

namespace NuTo
{
/**
 * @brief Engineering stress
 *
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
class EngineeringStress : public Eigen::Matrix<double, Voigt::Dim(TDim), 1>
{
    using Parent = Eigen::Matrix<double, Voigt::Dim(TDim), 1>;

public:
    EngineeringStress() = default;

    //! @brief This constructor allows you to construct EngineeringStress from Eigen expressions
    template <typename OtherDerived>
    EngineeringStress(const Eigen::MatrixBase<OtherDerived>& other)
        : Parent(other)
    {
    }

    //! @brief This constructor allows you to construct EngineeringStress from Eigen expressions
    template <typename OtherDerived>
    EngineeringStress(Eigen::MatrixBase<OtherDerived>&& other)
        : Parent(other)
    {
    }

    //! @brief This method allows you to assign Eigen expressions to EngineeringStress
    template <typename OtherDerived>
    EngineeringStress& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Parent::operator=(other);
        return *this;
    }

    //! @brief This method allows you to assign Eigen expressions to EngineeringStress
    template <typename OtherDerived>
    EngineeringStress& operator=(Eigen::MatrixBase<OtherDerived>&& other)
    {
        this->Parent::operator=(other);
        return *this;
    }
};
} /* namespace NuTo */
