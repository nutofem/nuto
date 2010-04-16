// $Id$

#ifndef DEFORMATIONGRADIENT2D_H_
#define DEFORMATIONGRADIENT2D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class DeformationGradient1D;
class DeformationGradient3D;
class EngineeringStrain2D;
class EngineeringStrain3D;
class GreenLagrangeStrain2D;
class GreenLagrangeStrain3D;
class Plane;

//! @brief ... two-dimensional deformation gradient
/*!
    In the two-dimensional case, deformation gradient \f$\boldsymbol{F}\f$ is defined as
    \f[\boldsymbol{F} = \begin{bmatrix}
      \dfrac{\partial x_1'}{\partial x_1} & \dfrac{\partial x_1'}{\partial x_2} \\
      \dfrac{\partial x_2'}{\partial x_1} & \dfrac{\partial x_2'}{\partial x_2}
    \end{bmatrix}, \f]
    where \f$\boldsymbol{x'}\f$ are the cartesian coordinates of the material point in the deformed configuration and
    \f$\boldsymbol{x}\f$ are the cartesian coordinates of the material point in the undeformed/initial configuration.
*/
//! @author Stefan Eckardt, ISM
//! @date November 2009
class DeformationGradient2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::Plane;
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    DeformationGradient2D();

    //! @brief ... copy constructor
    DeformationGradient2D(const DeformationGradient1D& rOther);

    //! @brief ... copy constructor
    DeformationGradient2D(const DeformationGradient2D& rOther);

    //! @brief ... get number of strain components
    //! @return ... number of strain components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get deformation gradient
    //! @return ... two-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    const double* GetDeformationGradient2D() const;

    //! @brief ... get deformation gradient
    //! @param  rDeformationGradient ... two-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    void GetDeformationGradient(NuTo::DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... get the three-dimensional deformation gradient
    /*!
     *  The two-dimensional deformation gradient is mapped to the three-dimensional case
     *  \f{align*}{
     *    \boldsymbol{F} = \begin{bmatrix}
      \dfrac{\partial x_1'}{\partial x_1} & \dfrac{\partial x_1'}{\partial x_2} & 0\\
      \dfrac{\partial x_2'}{\partial x_1} & \dfrac{\partial x_2'}{\partial x_2} & 0\\
      0 & 0 & 0
    \end{bmatrix}.
     *  \f}
     */
    //! @param  rDeformationGradient ... three-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    void GetDeformationGradient(NuTo::DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... set deformation gradient
    //! @param rDeformationGradient ... deformation gradient (column major format)
    //! @sa mDeformationGradient
    void SetDeformationGradient2D(const double* rDeformationGradient);

    //! @brief ... calculate reduced two-dimensional engineering (Cauchy) strain from deformation gradient
    /*!
     *  In the two-dimensional case the engineering strain reduces to
     *  \f{align*}{
     *      \varepsilon_{xx} &= \dfrac{\partial x_1'}{\partial x_1} - 1\\
     *      \varepsilon_{yy} &= \dfrac{\partial x_2'}{\partial x_2} - 1\\
     *      \varepsilon_{xy} &= \varepsilon_{yx} = \dfrac{1}{2} \left[ \dfrac{\partial x_1'}{\partial x_2} + \dfrac{\partial x_2'}{\partial x_1}\right].
     *  \f}
     *  Due to the symmetry of the second order engineering strain tensor the engineering strain components can be stored as vector
     *  \f[
     *       \boldsymbol{\varepsilon} = \begin{bmatrix}
     *         \varepsilon_{xx}\\
     *         \varepsilon_{yy}\\
     *         2 \varepsilon_{xy} \end{bmatrix} = \begin{bmatrix}
     *         \varepsilon_{xx}\\
     *         \varepsilon_{yy}\\
     *         \gamma_{xy} \end{bmatrix}.
     *  \f]
     */
    //! @param rEngineeringStrain ... engineering strain components
    void GetEngineeringStrain(NuTo::EngineeringStrain2D& rEngineeringStrain) const;

    //! @brief ... calculate three-dimensional engineering (Cauchy) strain from two-dimensional deformation gradient
    /*!
     *  The two-dimensional strain tensor is mapped to the three-dimensional case. The mapped engineering strain tensor reads
     *  \f{align*}{
     *      \boldsymbol{\varepsilon} = \begin{bmatrix}
     *        \varepsilon_{xx} & \varepsilon_{xy} & 0\\
     *        \varepsilon_{xy} & \varepsilon_{yy} & 0\\
     *        0 & 0 & 0
     *      \end{bmatrix}.
     *  \f}
     *  Due to the symmetry of the second order engineering strain tensor the engineering strain components can be stored as vector
     *  \f[
     *       \boldsymbol{\varepsilon} = \begin{bmatrix}
     *         \varepsilon_{xx}\\
     *         \varepsilon_{yy}\\
     *         0\\
     *         2 \varepsilon_{xy}\\
     *         0\\
     *         0
     *         \end{bmatrix} = \begin{bmatrix}
     *         \varepsilon_{xx}\\
     *         \varepsilon_{yy}\\
     *         0\\
     *         \gamma_{xy}\\
     *         0\\
     *         0
     *         \end{bmatrix}.
     *  \f]
     */
    //! @param rEngineeringStrain ... engineering strain components
    void GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const;

    //! @brief ... calculate reduced two-dimensional Green (Lagrangian) strain from deformation gradient
    /*!
     *  In the two-dimensional case, the Green (Lagrangian) strain reduces to
     *  \f{align*}{
     *      E_{xx} &= \dfrac{1}{2}\left[\left(\dfrac{\partial x_1'}{\partial x_1}\right)^2 +  \left(\dfrac{\partial x_2'}{\partial x_1}\right)^2 - 1\right]\\
     *      E_{yy} &= \dfrac{1}{2}\left[\left(\dfrac{\partial x_1'}{\partial x_2}\right)^2 +  \left(\dfrac{\partial x_2'}{\partial x_2}\right)^2 - 1\right]\\
     *      E_{xy} &= E_{yx} = \dfrac{1}{2} \left[\dfrac{\partial x_1'}{\partial x_1} \dfrac{\partial x_1'}{\partial x_2} + \dfrac{\partial x_2'}{\partial x_2} \dfrac{\partial x_2'}{\partial x_1} \right].
     *  \f}
     *  Due to the symmetry of the second order Green strain tensor the Green strain components can be stored as vector
     *  \f[
     *      \boldsymbol{E} = \begin{bmatrix}
     *        E_{xx}\\
     *        E_{yy}\\
     *        2 E_{xy}
     *      \end{bmatrix}.
     *  \f]
     */
    //! @param rGreenStrain ... Green strain components
    void GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain2D& rGreenLagrangeStrain) const;

    //! @brief ... calculate three-dimensional Green (Lagrangian) strain from two-dimensional deformation gradient
    /*!
     *  The two-dimensional Green strain tensor is mapped to the three-dimensional case. The mapped Green strain tensor reads
     *  \f{align*}{
     *      \boldsymbol{E} = \begin{bmatrix}
     *        E_{xx} & E_{xy} & 0\\
     *        E_{xy} & E_{yy} & 0\\
     *        0 & 0 & 0
     *      \end{bmatrix}.
     *  \f}
     *  Due to the symmetry of the second order Green strain tensor the Green strain components can be stored as vector
     *  \f[
     *       \boldsymbol{E} = \begin{bmatrix}
     *         E_{xx}\\
     *         E_{yy}\\
     *         0\\
     *         2 E_{xy}\\
     *         0\\
     *         0
     *         \end{bmatrix}.
     *  \f]
     */
    //! @param rGreenLagrangeStrain ... green strain components
    void GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D& rGreenStrain) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(this->mDeformationGradient);
    }
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... deformation gradient
    //!
    //! The deformation gradient is stored in column major format:
    //! \f$\left[\dfrac{\partial x_1'}{\partial x_1},\dfrac{\partial x_2'}{\partial x_1},\dfrac{\partial x_1'}{\partial x_2},\dfrac{\partial x_2'}{\partial x_2}\right]\f$
    //!
    double mDeformationGradient[4];
};

}

#endif // DEFORMATIONGRADIENT2D_H_
