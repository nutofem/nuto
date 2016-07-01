#pragma once

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include <type_traits>

namespace NuTo
{

//! @brief ... Engineering strain
/*!
 *  Due to the symmetry of the second order engineering strain tensor the engineering strain components can be stored as vector
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
class EngineeringStrain : public ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)>
{
public:
    EngineeringStrain()                                       = default;
    EngineeringStrain(const EngineeringStrain& )              = default;
    EngineeringStrain(      EngineeringStrain&&)              = default;

    EngineeringStrain& operator=(const EngineeringStrain& )   = default;
    EngineeringStrain& operator=(      EngineeringStrain&&)   = default;

    virtual std::unique_ptr<ConstitutiveIOBase> clone() override;
    //! @brief returns I1 - the first strain invariant of the characteristic equation
    //! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
    //! @remark only implemented for TDim = 3
    template<int U = TDim, typename std::enable_if<U == 3, int>::type = 0>
    double InvariantI1() const;

    //! @brief returns I2 - the first strain invariant of the characteristic equation
    //! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
    //! @remark only implemented for TDim = 3
    template<int U = TDim, typename std::enable_if<U == 3, int>::type = 0>
    double InvariantI2() const;

    //! @brief returns I3 - the first strain invariant of the characteristic equation
    //! \f[ \lambda^3 - I_1 \lambda^2 + I_2 \lambda - I_3  \f]
    //! @remark only implemented for TDim = 3
    template<int U = TDim, typename std::enable_if<U == 3, int>::type = 0>
    double InvariantI3() const;

    //! @brief returns J2 - the second deviatoric strain invariant of the characteristic equation
    //! \f[ \lambda^3 - J_1 \lambda^2 - J_2 \lambda - J_3  \f] Note the minus sign in front of J2. This is not consistent with the I2 invariant but apparently common practice.
    //! @remark only implemented for TDim = 3
    template<int U = TDim, typename std::enable_if<U == 3, int>::type = 0>
    double InvariantJ2() const;


    //! @brief calculates the 3D engineering strain.
    //! @param rNu ... poisson ratio
    //! @param rSectionType ... defaulted to dummy, needed in 2D for plane stress / plane strain
    EngineeringStrain<3> As3D(double rNu, Section::eSectionType rSectionType = Section::VOLUME) const;

    //! @brief returns the deviatoric part
    //! @remark only implemented for TDim = 3
    template<int U = TDim, typename std::enable_if<U == 3, int>::type = 0>
    EngineeringStrain<3> Deviatoric() const;

    const EngineeringStrain<1>& AsEngineeringStrain1D() const override {throw;}
    const EngineeringStrain<2>& AsEngineeringStrain2D() const override {throw;}
    const EngineeringStrain<3>& AsEngineeringStrain3D() const override {throw;}
    EngineeringStrain<1>& AsEngineeringStrain1D() override {throw;}
    EngineeringStrain<2>& AsEngineeringStrain2D() override {throw;}
    EngineeringStrain<3>& AsEngineeringStrain3D() override {throw;}

};

template<> inline const NuTo::EngineeringStrain<1>& NuTo::EngineeringStrain<1>::AsEngineeringStrain1D() const {return *this;}
template<> inline const NuTo::EngineeringStrain<2>& NuTo::EngineeringStrain<2>::AsEngineeringStrain2D() const {return *this;}
template<> inline const NuTo::EngineeringStrain<3>& NuTo::EngineeringStrain<3>::AsEngineeringStrain3D() const {return *this;}
template<> inline NuTo::EngineeringStrain<1>& NuTo::EngineeringStrain<1>::AsEngineeringStrain1D() {return *this;}
template<> inline NuTo::EngineeringStrain<2>& NuTo::EngineeringStrain<2>::AsEngineeringStrain2D() {return *this;}
template<> inline NuTo::EngineeringStrain<3>& NuTo::EngineeringStrain<3>::AsEngineeringStrain3D() {return *this;}

} /* namespace NuTo */


