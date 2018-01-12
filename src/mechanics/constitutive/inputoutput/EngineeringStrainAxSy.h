#pragma once

#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"

#include <type_traits>

namespace NuTo
{

//! @brief ... Engineering strain
/*!
 *  Due to the symmetry of the second order engineering strain tensor the engineering strain components can be stored as
 * vector
 *  \f[
 *       \boldsymbol{\varepsilon} = \begin{bmatrix}
 *         \varepsilon_{rr}\\
 *         \varepsilon_{zz}\\
 *         \varepsilon_{thth}\\
 *         2 \varepsilon_{rz}
 *       \end{bmatrix} = \begin{bmatrix}
 *         \varepsilon_{r}\\
 *         \varepsilon_{z}\\
 *         \varepsilon_{th}\\
 *         \gamma_{rz}
 *       \end{bmatrix}.
 *  \f]
 *
 */
class EngineeringStrainAxSy : public ConstitutiveVector<4>
{
public:
    EngineeringStrainAxSy() = default;
    EngineeringStrainAxSy(std::initializer_list<double> initList);
    EngineeringStrainAxSy(const EngineeringStrainAxSy&) = default;
    EngineeringStrainAxSy(EngineeringStrainAxSy&&) = default;

    EngineeringStrainAxSy& operator=(const EngineeringStrainAxSy&) = default;
    EngineeringStrainAxSy& operator=(EngineeringStrainAxSy&&) = default;

    //! @brief used in the EngineeringStrain class, not used here
    virtual std::unique_ptr<ConstitutiveIOBase> clone() override
    {
    	return std::make_unique<EngineeringStrainAxSy>(*this);
    }

    //! @brief returns I1 - the first strain invariant of the characteristic equation
    //! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
    //! @remark only implemented for TDim = 3
    double InvariantI1() const;

    //! @brief returns I2 - the first strain invariant of the characteristic equation
    //! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
    //! @remark only implemented for TDim = 3
    double InvariantI2() const;

    //! @brief returns I3 - the first strain invariant of the characteristic equation
    //! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
    //! @remark only implemented for TDim = 3
    double InvariantI3() const;

    //! @brief returns J2 - the second deviatoric strain invariant of the characteristic equation
    //! \f[ \lambda^3 - J_1 \lambda^2 - J_2 \lambda - J_3  \f] Note the minus sign in front of J2. This is not
    //! consistent with the I2 invariant but apparently common practice.
    //! @remark only implemented for TDim = 3
    double InvariantJ2() const;


    //! @brief Calculates the 3D engineering strain.
    //! @param rNu Poisson ratio.
    //! @param rSectionType Defaults to plane stress; needed in 2D.
    EngineeringStrain<3> As3D() const;

    //! @brief returns the deviatoric part
    //! @remark only implemented for TDim = 3
    EngineeringStrainAxSy Deviatoric() const;

    const EngineeringStrainAxSy& AsEngineeringStrainAxSy() const override
    {
        return *this;
    }
    EngineeringStrainAxSy& AsEngineeringStrainAxSy() override
    {
        return *this;
    }
};

//inline const NuTo::EngineeringStrainAxSy& NuTo::EngineeringStrainAxSy::AsEngineeringStrainAxSy() const
//{
//    return *this;
//}
//inline NuTo::EngineeringStrainAxSy& NuTo::EngineeringStrainAxSy::AsEngineeringStrainAxSy()
//{
//    return *this;
//}

} /* namespace NuTo */
