#pragma once

#include <eigen3/Eigen/Core>
#include "nuto/mechanics/constitutive/Voigt.h"

namespace NuTo
{
/**
 * Engineering strain
 *
 *  Due to the symmetry of the second order engineering strain tensor the engineering strain components can be stored as
 * vector
 *  \f[
 *       \boldsymbol{\varepsilon} = \begin{bmatrix}
 *         \varepsilon_{xx}\\
 *         \varepsilon_{yy}\\
 *         \varepsilon_{zz}\\
 *         2 \varepsilon_{xy}\\
 *         2 \varepsilon_{yz}\\
 *         2 \varepsilon_{zx}
 *       \end{bmatrix} = \begin{bmatrix}
 *         \varepsilon_{x}\\
 *         \varepsilon_{y}\\
 *         \varepsilon_{z}\\
 *         \gamma_{xy} \\
 *         \gamma_{yz} \\
 *         \gamma_{zx}
 *       \end{bmatrix}.
 *  \f]
 *
 */
template <int TDim>
class EngineeringStrain : public Eigen::Matrix<double, Voigt::Dim(TDim), 1>
{
    using Parent = Eigen::Matrix<double, Voigt::Dim(TDim), 1>;

public:
    EngineeringStrain() = default;

    //! @brief This constructor allows you to construct EngineeringStrain from Eigen expressions
    template <typename OtherDerived>
    EngineeringStrain(const Eigen::MatrixBase<OtherDerived>& other)
        : Parent(other)
    {
    }

    //! @brief This constructor allows you to construct EngineeringStrain from Eigen expressions
    template <typename OtherDerived>
    EngineeringStrain(Eigen::MatrixBase<OtherDerived>&& other)
        : Parent(other)
    {
    }

    //! @brief This method allows you to assign Eigen expressions to EngineeringStrain
    template <typename OtherDerived>
    EngineeringStrain& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Parent::operator=(other);
        return *this;
    }

    //! @brief This method allows you to assign Eigen expressions to EngineeringStrain
    template <typename OtherDerived>
    EngineeringStrain& operator=(Eigen::MatrixBase<OtherDerived>&& other)
    {
        this->Parent::operator=(other);
        return *this;
    }
};
} /* namespace NuTo */
