// $Id$

#ifndef DEFORMATIONGRADIENT3D_H_
#define DEFORMATIONGRADIENT3D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class DeformationGradient1D;
class DeformationGradient2D;
class EngineeringStrain3D;
class GreenLagrangeStrain3D;
class Solid;
//! @brief ... three-dimensional deformation gradient
/*!
   In the two-dimensional case, the deformation gradient \f$\boldsymbol{F}\f$ is defined as
   \f[\boldsymbol{F} = \begin{bmatrix}
       \dfrac{\partial x_1'}{\partial x_1} & \dfrac{\partial x_1'}{\partial x_2} & \dfrac{\partial x_1'}{\partial x_3} \\
       \dfrac{\partial x_2'}{\partial x_1} & \dfrac{\partial x_2'}{\partial x_2} & \dfrac{\partial x_2'}{\partial x_3} \\
       \dfrac{\partial x_3'}{\partial x_1} & \dfrac{\partial x_3'}{\partial x_2} & \dfrac{\partial x_3'}{\partial x_3}
      \end{bmatrix}, \f]
    where \f$\boldsymbol{x'}\f$ are the cartesian coordinates of the material point in the deformed configuration and
    \f$\boldsymbol{x}\f$ are the cartesian coordinates of the material point in the undeformed/initial configuration.
*/
//! @author Stefan Eckardt, ISM
//! @date November 2009
class DeformationGradient3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::Solid;
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    DeformationGradient3D();

    //! @brief ... copy constructor
    DeformationGradient3D(const DeformationGradient1D& rOther);

    //! @brief ... copy constructor
    DeformationGradient3D(const DeformationGradient2D& rOther);

    //! @brief ... copy constructor
    DeformationGradient3D(const DeformationGradient3D& rOther);

    //! @brief ... get number of strain components
    //! @return ... number of strain components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get deformation gradient
    //! @return ... deformation gradient (column major format)
    //! @sa mDeformationGradient
    const double* GetDeformationGradient3D() const;

    //! @brief ... get deformation gradient
    //! @param rDeformationGradient ... deformation gradient (column major format)
    //! @sa mDeformationGradient
    void GetDeformationGradient(NuTo::DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... set deformation gradient
    //! @param rDeformationGradient ... deformation gradient (column major format)
    //! @sa mDeformationGradient
    void SetDeformationGradient3D(const double* rDeformationGradient);

    //! @brief ... calculate engineering (Cauchy) strain from deformation gradient
    /*!
     *  In the three-dimensional case the components of the engineering strain tensor read
     *  \f{align*}{
     *      \varepsilon_{xx} &= \dfrac{\partial x_1'}{\partial x_1} - 1\\
     *      \varepsilon_{yy} &= \dfrac{\partial x_2'}{\partial x_2} - 1\\
     *      \varepsilon_{zz} &= \dfrac{\partial x_3'}{\partial x_3} - 1\\
     *      \varepsilon_{xy} &= \varepsilon_{yx} = \dfrac{1}{2} \left[ \dfrac{\partial x_1'}{\partial x_2} + \dfrac{\partial x_2'}{\partial x_1}\right]\\
     *      \varepsilon_{yz} &= \varepsilon_{zy} = \dfrac{1}{2} \left[ \dfrac{\partial x_2'}{\partial x_3} + \dfrac{\partial x_3'}{\partial x_2}\right]\\
     *      \varepsilon_{zx} &= \varepsilon_{xz} = \dfrac{1}{2} \left[ \dfrac{\partial x_1'}{\partial x_3} + \dfrac{\partial x_3'}{\partial x_1}\right].
     *  \f}
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
     *         \varepsilon_{xx}\\
     *         \varepsilon_{yy}\\
     *         \varepsilon_{zz}\\
     *         \gamma_{xy} \\
     *         \gamma_{yz} \\
     *         \gamma_{zx}
     *       \end{bmatrix}.
     *  \f]
     */
    //! @param rEngineeringStrain ... engineering strain components
    void GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const;

    //! @brief ... calculate Green (Lagrangian) strain from deformation gradient
    /*!
     *  In the three-dimensional case, the components of the Green (Lagrangian) strain tensor read
     *  \f{align*}{
     *      E_{xx} &= \dfrac{1}{2}\left[\left(\dfrac{\partial x_1'}{\partial x_1}\right)^2 +  \left(\dfrac{\partial x_2'}{\partial x_1}\right)^2 +  \left(\dfrac{\partial x_3'}{\partial x_1}\right)^2- 1\right]\\
     *      E_{yy} &= \dfrac{1}{2}\left[\left(\dfrac{\partial x_1'}{\partial x_2}\right)^2 +  \left(\dfrac{\partial x_2'}{\partial x_2}\right)^2 +  \left(\dfrac{\partial x_3'}{\partial x_2}\right)^2- 1\right]\\
     *      E_{zz} &= \dfrac{1}{2}\left[\left(\dfrac{\partial x_1'}{\partial x_3}\right)^2 +  \left(\dfrac{\partial x_2'}{\partial x_3}\right)^2 +  \left(\dfrac{\partial x_3'}{\partial x_3}\right)^2- 1\right]\\
     *      E_{xy} &= E_{yx} = \dfrac{1}{2} \left[\dfrac{\partial x_1'}{\partial x_1} \dfrac{\partial x_1'}{\partial x_2} + \dfrac{\partial x_2'}{\partial x_1} \dfrac{\partial x_2'}{\partial x_2} +  \dfrac{\partial x_3'}{\partial x_2} \dfrac{\partial x_3'}{\partial x_1}\right]\\
     *      E_{yz} &= E_{zy} = \dfrac{1}{2} \left[\dfrac{\partial x_1'}{\partial x_3} \dfrac{\partial x_1'}{\partial x_2} + \dfrac{\partial x_2'}{\partial x_3} \dfrac{\partial x_2'}{\partial x_2} +  \dfrac{\partial x_3'}{\partial x_3} \dfrac{\partial x_3'}{\partial x_2}\right]\\
     *      E_{zx} &= E_{xz} = \dfrac{1}{2} \left[\dfrac{\partial x_1'}{\partial x_3} \dfrac{\partial x_1'}{\partial x_1} + \dfrac{\partial x_2'}{\partial x_3} \dfrac{\partial x_2'}{\partial x_1} +  \dfrac{\partial x_3'}{\partial x_3} \dfrac{\partial x_3'}{\partial x_1}\right].
     *  \f}
     *  Due to the symmetry of the second order Green strain tensor the Green strain components can be stored as vector
     *  \f[
     *      \boldsymbol{E} = \begin{bmatrix}
     *        E_{xx}\\
     *        E_{yy}\\
     *        E_{zz}\\
     *        2 E_{xy}\\
     *        2 E_{yz}\\
     *        2 E_{zx}
     *      \end{bmatrix}.
     *  \f]
     */
    //! @param rGreenStrain ... Green strain components
    void GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D& rGreenStrain) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... deformation gradient
    //!
    //! The deformation gradient is stored in column major format:
    //! \f$\left[\dfrac{\partial x_1'}{\partial x_1},\dfrac{\partial x_2'}{\partial x_1},\dfrac{\partial x_3'}{\partial x_1},\dfrac{\partial x_1'}{\partial x_2},\dfrac{\partial x_2'}{\partial x_2},\dfrac{\partial x_3'}{\partial x_2},\dfrac{\partial x_1'}{\partial x_3},\dfrac{\partial x_2'}{\partial x_3},\dfrac{\partial x_3'}{\partial x_3}\right]\f$
    //!
    double mDeformationGradient[9];
};

}

#endif // DEFORMATIONGRADIENT3D_H_
