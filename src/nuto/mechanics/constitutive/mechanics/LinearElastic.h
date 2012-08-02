// $Id: LinearElastic.h 112 2009-11-17 16:43:15Z unger3 $

#ifndef LINEARELASTIC_H
#define LINEARELASTIC_H

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutivePiolaKirchhoffIIGreenLagrange.h"


namespace NuTo
{
class ConstitutiveInput3D;
class ConstitutiveOutput3D;
//! @brief ... linear elastic material model
/*!
 * Assuming linear elastic material behavior, the one-dimensional constitutive relationship reads
 * \f{align*}{
 *   \sigma_{xx} = E \varepsilon_{xx}
 * \f}
 * where \f$ E \f$ is the Young's modulus, \f$ \sigma_{xx} \f$ is the Cauchy stress
 * and \f$ \varepsilon_{xx} \f$ is the Engineering strain vector.
 * It is to be noted, that this relationship holds for Green strains and the second Piola-Kirchhoff stress as well.
 * Assuming linear elastic isotropic material behavior, the general three-dimensional constitutive relationship reads
 * \f{align*}{
 *  \begin{bmatrix}
 *    \sigma_{xx}\\
 *    \sigma_{yy}\\
 *    \sigma_{zz}\\
 *    \sigma_{xy}\\
 *    \sigma_{yz}\\
 *    \sigma_{zx}
 *  \end{bmatrix} = \dfrac{E (1 - \nu)}{(1+\nu)(1-2\nu)} \begin{bmatrix}
 *    1 & \dfrac{\nu}{1-\nu} & \dfrac{\nu}{1-\nu} & 0 & 0 & 0\\
 *    \dfrac{\nu}{1-\nu} & 1 & \dfrac{\nu}{1-\nu} & 0 & 0 & 0\\
 *    \dfrac{\nu}{1-\nu} & \dfrac{\nu}{1-\nu} & 1 & 0 & 0 & 0\\
 *    0 & 0 & 0 & \dfrac{1-2\nu}{2(1-\nu)} & 0 & 0\\
 *    0 & 0 & 0 & 0 & \dfrac{1-2\nu}{2(1-\nu)} & 0\\
 *    0 & 0 & 0 & 0 & 0 & \dfrac{1-2\nu}{2(1-\nu)}
 *  \end{bmatrix}
 *  \begin{bmatrix}
 *    \varepsilon_{xx}\\
 *    \varepsilon_{yy}\\
 *    \varepsilon_{zz}\\
 *    \gamma_{xy}\\
 *    \gamma_{yz}\\
 *    \gamma_{zx}
 *  \end{bmatrix},
 * \f}
 * two-dimensional plain stress constitutive relationship reads
 * \f{align*}{
 *  \begin{bmatrix}
 *    \sigma_{xx}\\
 *    \sigma_{yy}\\
 *    \sigma_{xy}
 *  \end{bmatrix} = \dfrac{E}{(1+\nu^2)} \begin{bmatrix}
 *    1 & \nu & 0\\
 *    \nu & 1 & 0\\
 *    0 & 0 \dfrac{1-\nu}{2}
 *  \end{bmatrix}
 *  \begin{bmatrix}
 *    \varepsilon_{xx}\\
 *    \varepsilon_{yy}\\
 *    \gamma_{xy}
 *  \end{bmatrix},
 * \f} 
 * where \f$ E \f$ is the Young's modulus, \f$ \nu \f$ is the Poisson's ratio,
 * \f$ \boldsymbol{\sigma} \f$ are the components of the Cauchy stress vector,
 * and \f$ \boldsymbol{\varepsilon} \f$ are the components of the Engineering strain vector.
 * It is to be noted, that this relationship holds for Green strains and the second Piola-Kirchhoff stress as well.
 */
//! @author Jörg F. Unger, ISM
//! @date November 2009
class LinearElastic: public ConstitutiveEngineeringStressStrain, public ConstitutivePiolaKirchhoffIIGreenLagrange
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    LinearElastic();

	// Engineering strain /////////////////////////////////////
	//! @brief ... calculate 3D engineering strain vector from 2D engineering strain vector
	//! @param rElement ... element
	//! @param rIp ... integration point
	//! @param rEngineeringStrain2D ... engineering strain
	//! @param rEngineeringStrain3D ... engineering strain
	Error::eError GetEngineeringStrainFromEngineeringStrain(const ElementBase* rElement, int rIp,
				  const EngineeringStrain2D& rEngineeringStrain2D, EngineeringStrain3D& rEngineeringStrain3D) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    Error::eError GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    Error::eError GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    Error::eError GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;


    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rCauchyStress ... Cauchy stress
    Error::eError GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, double& rDamage) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, double& rDamage) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, double& rDamage) const;


    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetInternalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const;

    //! @brief ... calculates the difference of the elastic strain between the current state and the previous update
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
    void GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const;
    ///////////////////////////////////////////////////////////////////////////

    // Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
    //! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
    Error::eError GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, SecondPiolaKirchhoffStress1D& rSecondPiolaKirchhoffStress) const;

    // Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
    //! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
    Error::eError GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, SecondPiolaKirchhoffStress2D& rSecondPiolaKirchhoffStress) const;

    // Second PiolaKirchhoff Stress - Green Lagrange Strains /////////////////////////////////////
    //! @brief ... calculate second Piola-Kirchhoff stress from GreenLagrangeStrains (which is calculated from the deformation gradient)
   //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rSecondPiolaKirchhoffStress ... Second Piola Kirchhoff Stress
    Error::eError GetSecondPiolaKirchhoffStressFromGreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient3D& rDeformationGradient, SecondPiolaKirchhoffStress3D& rSecondPiolaKirchhoffStress) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient,
            ConstitutiveTangentLocal1x1& rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient,
            ConstitutiveTangentLocal3x3& rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient,
            ConstitutiveTangentLocal6x6& rTangent) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError UpdateStaticData_SecondPiolaKirchhoffStress_GreenLagrangeStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain1D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain2D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataSecondPiolaKirchhoffStress_GreenLagrangeStrain3D(const ElementBase* rElement) const;


    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetInternalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetInternalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetInternalEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    Error::eError GetElasticEnergy_SecondPiolaKirchhoffStress_GreenLagrangeStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient, double& rEnergy) const;
    ///////////////////////////////////////////////////////////////////////////

    // calculate coefficients of the material matrix
    void CalculateCoefficients2DPlainStress(double& C11, double& C12, double& C33) const;
    void CalculateCoefficients3D(double& C11, double& C12, double& C44) const;

    // parameters /////////////////////////////////////////////////////////////
    //! @brief ... get density
    //! @return ... density
    virtual double GetDensity() const;

    //! @brief ... set density
    //! @param rRho ... density
    virtual void SetDensity(double rRho);

    //! @brief ... get Young's modulus
    //! @return ... Young's modulus
    double GetYoungsModulus() const;

    //! @brief ... set Young's modulus
    //! @param rE ... Young's modulus
    void SetYoungsModulus(double rE);

    //! @brief ... get Poisson's ratio
    //! @return ... Poisson's ratio
    double GetPoissonsRatio() const;

    //! @brief ... set Poisson's ratio
    //! @param rNu ... Poisson's ratio
    void SetPoissonsRatio(double rNu);

    //! @brief ... get thermal expansion coefficient
    //! @return ... thermal expansion coefficient
    double GetThermalExpansionCoefficient() const;

    //! @brief ... set thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void SetThermalExpansionCoefficient(double rNu);
    ///////////////////////////////////////////////////////////////////////////

    //! @brief ... get dimension of the constitutive relationship
    //! @return ... dimension of the constitutive relationship (1, 2 or 3)
    int GetGlobalDimension() const;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters()const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    //! @param rLogger stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const
    {
    	return false;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... density \f$ \rho \f$
    double mRho;

    //! @brief ... thermal expansion coefficient \f$ \alpha \f$
    double mThermalExpansionCoefficient;

    //! @brief ... check if density is positive
    //! @param rRho ... density
    void CheckDensity(double rRho) const;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void CheckThermalExpansionCoefficient(double rAlpha) const;
};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LinearElastic)

//this is due to the diamond structure (virtual public derivative)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::ConstitutiveEngineeringStressStrain, NuTo::LinearElastic>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::ConstitutivePiolaKirchhoffIIGreenLagrange, NuTo::LinearElastic>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif // LINEARELASTIC_H_
