// $Id$

#ifndef DEFORMATIONGRADIENT1D_H_
#define DEFORMATIONGRADIENT1D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class DeformationGradient2D;
class DeformationGradient3D;
class EngineeringStrain1D;
class EngineeringStrain2D;
class EngineeringStrain3D;
class GreenLagrangeStrain1D;
class GreenLagrangeStrain2D;
class GreenLagrangeStrain3D;
class Truss;

//! @brief ... one-dimensional deformation gradient
/*!
    In the one-dimensional case, deformation gradient \f$\boldsymbol{F}\f$ is defined as
    \f[\boldsymbol{F} = \begin{bmatrix}
      \dfrac{\partial x_1'}{\partial x_1}
    \end{bmatrix}, \f]
    where \f$\boldsymbol{x'}\f$ are the cartesian coordinates of the material point in the deformed configuration and
    \f$\boldsymbol{x}\f$ are the cartesian coordinates of the material point in the undeformed/initial configuration.
*/
//! @author Stefan Eckardt, ISM
//! @date November 2009
class DeformationGradient1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::Truss;
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    DeformationGradient1D();

    //! @brief ... copy constructor
    DeformationGradient1D(const DeformationGradient1D& rOther);

    //! @brief ... get number of strain components
    //! @return ... number of strain components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get deformation gradient
    //! @return ... one-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    const double* GetDeformationGradient1D() const;

    //! @brief ... get deformation gradient
    //! @param  rDeformationGradient ... one-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    void GetDeformationGradient(NuTo::DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... get deformation gradient
    /*!
     *  The mapped two-dimensional deformation gradient reads
     *  \f{align*}{
     *    \boldsymbol{F} = \begin{bmatrix}
      \dfrac{\partial x_1'}{\partial x_1} & 0\\
      0 & 0
    \end{bmatrix}.
     *  \f}
     */
    //! @param  rDeformationGradient ... two-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    void GetDeformationGradient(NuTo::DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... get deformation gradient
    /*!
     *  The mapped two-dimensional deformation gradient reads
     *  \f{align*}{
     *    \boldsymbol{F} = \begin{bmatrix}
      \dfrac{\partial x_1'}{\partial x_1} & 0 & 0\\
      0 & 0 & 0\\
      0 & 0 & 0
    \end{bmatrix}.
     *  \f}
     */
    //! @param  rDeformationGradient ... three-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    void GetDeformationGradient(NuTo::DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... set deformation gradient
    //! @param rDeformationGradient ... one-dimensional deformation gradient (column major format)
    //! @sa mDeformationGradient
    void SetDeformationGradient1D(const double* rDeformationGradient);

    //! @brief ... calculate reduced one-dimensional engineering (Cauchy) strain from deformation gradient
    /*!
     *  In the one-dimensional case the engineering strain reduces to
     *  \f[ \varepsilon_{xx} = \dfrac{\partial x_1'}{\partial x_1} - 1. \f]
     *  Due to the symmetry of the second order engineering strain tensor the engineering strain components can be stored as vector
     *  \f[ \boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_{xx} \end{bmatrix}. \f]
     */
    //! @param rEngineeringStrain ... engineering strain components
    void GetEngineeringStrain(NuTo::EngineeringStrain1D& rEngineeringStrain) const;

    //! @brief ... calculate reduced two-dimensional engineering (Cauchy) strain from one-dimensional deformation gradient
    /*!
     *  The one-dimensional engineering strain is mapped to the two-dimensional case.
     *  The mapped engineering strain tensor reads
     *  \f{align*}{
     *    \boldsymbol{\varepsilon} = \begin{bmatrix}
     *      \varepsilon_{xx} & 0\\
     *      0 & 0\\
     *    \end{bmatrix}.
     *  \f}
     *  Due to the symmetry of the engineering strain tensor the engineering strain components are stored as vector
     *  \f{align*}{
     *    \boldsymbol{\varepsilon} = \begin{bmatrix}
     *      \varepsilon_{xx}\\
     *      0\\
     *      0
     *    \end{bmatrix}.
     *  \f}
     */
    //! @param rEngineeringStrain ... engineering strain components
    void GetEngineeringStrain(NuTo::EngineeringStrain2D& rEngineeringStrain) const;

    //! @brief ... calculate three-dimensional engineering (Cauchy) strain from one-dimensional deformation gradient
    /*!
     *  The one-dimensional engineering strain is mapped to the three-dimensional case.
     *  The mapped engineering strain tensor reads
     *  \f{align*}{
     *    \boldsymbol{\varepsilon} = \begin{bmatrix}
     *      \varepsilon_{xx} & 0 & 0\\
     *      0 & 0 & 0\\
     *      0 & 0 & 0
     *    \end{bmatrix}.
     *  \f}
     *  Due to the symmetry of the engineering strain tensor the engineering strain components are stored as vector
     *  \f{align*}{
     *    \boldsymbol{\varepsilon} = \begin{bmatrix}
     *      \varepsilon_{xx}\\
     *      0\\
     *      0\\
     *      0\\
     *      0\\
     *      0
     *    \end{bmatrix}.
     *  \f}
     */
    //! @param rEngineeringStrain ... engineering strain components
    void GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const;

    //! @brief ... calculate reduced one-dimensional Green (Lagrangian) strain from deformation gradient
    /*!
     *  In the one-dimensional case, the Green (Lagrangian) strain reduces to
     *  \f[ E_{xx} = \dfrac{1}{2}\left[\left(\dfrac{\partial x_1'}{\partial x_1}\right)^2  - 1\right]. \f]
     *  Due to the symmetry of the second order Green strain tensor the Green strain components can be stored as vector
     *  \f[ \boldsymbol{E} = \begin{bmatrix} E_{xx} \end{bmatrix}. \f]
     */
    //! @param rGreenLagrangeStrain ... Green strain components
    void GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain1D& rGreenLagrangeStrain) const;

    //! @brief ... calculate reduced two-dimensional Green (Lagrangian) strain from one-dimensional deformation gradient
    /*!
     *  The one-dimensional Green strain is mapped to the two-dimensional case.
     *  The mapped Green strain tensor reads
     *  \f{align*}{
     *    \boldsymbol{E} = \begin{bmatrix}
     *      E_{xx} & 0\\
     *      0 & 0\\
     *    \end{bmatrix}.
     *  \f}
     *  Due to the symmetry of the Green strain tensor the Green strain components are stored as vector
     *  \f{align*}{
     *    \boldsymbol{E} = \begin{bmatrix}
     *      E_{xx}\\
     *      0\\
     *      0
     *    \end{bmatrix}.
     *  \f}
     */
    //! @param rGreenStrain ... Green strain components
    void GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain2D&  rGreenLagrangeStrain) const;

    //! @brief ... calculate three-dimensional Green (Lagrangian) strain from one-dimensional deformation gradient
    /*!
     *  The one-dimensional Green strain is mapped to the three-dimensional case.
     *  The mapped Green strain tensor reads
     *  \f{align*}{
     *    \boldsymbol{E} = \begin{bmatrix}
     *      E_{xx} & 0 & 0\\
     *      0 & 0 & 0\\
     *      0 & 0 & 0
     *    \end{bmatrix}.
     *  \f}
     *  Due to the symmetry of the Green strain tensor the Green strain components are stored as vector
     *  \f{align*}{
     *    \boldsymbol{E} = \begin{bmatrix}
     *      E_{xx}\\
     *      0\\
     *      0\\
     *      0\\
     *      0\\
     *      0
     *    \end{bmatrix}.
     *  \f}
     */
    //! @param rGreenStrain ... Green strain components
    void GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D&  rGreenLagrangeStrain) const;

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
    //! The deformation gradient is stored in column major format: \f$\left[\dfrac{\partial x_1'}{\partial x_1}\right]\f$
    //!
    double mDeformationGradient;
};

}

#endif // DEFORMATIONGRADIENT1D_H_ 
