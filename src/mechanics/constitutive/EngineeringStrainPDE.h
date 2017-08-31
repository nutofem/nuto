#pragma once

#include <eigen3/Eigen/Core>
#include "mechanics/constitutive/Voigt.h"

namespace NuTo
{

//! @brief ... Engineering strain
/*!
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
class EngineeringStrainPDE : public Eigen::Matrix<double, Voigt::Dim(TDim), 1>
{
    using Parent = Eigen::Matrix<double, Voigt::Dim(TDim), 1>;
public:
    EngineeringStrainPDE() = default;

    // This constructor allows you to construct EngineeringStrainPDE from Eigen expressions
    template <typename OtherDerived>
    EngineeringStrainPDE(const Eigen::MatrixBase<OtherDerived>& other)
        : Parent(other) 
    {
    }
    
    // This constructor allows you to construct EngineeringStrainPDE from Eigen expressions
    template <typename OtherDerived>
    EngineeringStrainPDE(Eigen::MatrixBase<OtherDerived>&& other)
        : Parent(other) 
    {
    }
    
    // This method allows you to assign Eigen expressions to EngineeringStrainPDE
    template <typename OtherDerived>
    EngineeringStrainPDE& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Parent::operator=(other);
        return *this;
    }
    
    // This method allows you to assign Eigen expressions to EngineeringStrainPDE
    template <typename OtherDerived>
    EngineeringStrainPDE& operator=(Eigen::MatrixBase<OtherDerived>&& other)
    {
        this->Parent::operator=(other);
        return *this;
    }
};
} /* namespace NuTo */
