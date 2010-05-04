// $ld: $ 
#ifndef CONSTITUTIVENONLOCALDAMAGEPLASTICITY_H_
#define CONSTITUTIVENONLOCALDAMAGEPLASTICITY_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutivePiolaKirchhoffIIGreenLagrange.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 26, 2010
//! @brief ...
class NonlocalDamagePlasticity : public ConstitutiveEngineeringStressStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    NonlocalDamagePlasticity();

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering strain from deformation gradient in 1D (truss is assumed to be plane stress)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering strain from deformation gradient in 2D (e.g. plane stress)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering strain from deformation gradient in 3D
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;

    //  Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStrain ... engineering strain
    void GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const;


    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rEngineeringStress ... engineering stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    // Engineering stress - Engineering strain /////////////////////////////////////
    //! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rCauchyStress ... Cauchy stress
    void GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
    		      const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient,
            ConstitutiveTangentLocal1x1& rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient,
            ConstitutiveTangentLocal3x3& rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient,
            ConstitutiveTangentLocal6x6& rTangent) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
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
    double GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    double GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
            const DeformationGradient3D& rDeformationGradient) const;

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

    // calculate coefficients of the material matrix
    void CalculateCoefficients3D(double& C11, double& C12, double& C44) const;

    // parameters /////////////////////////////////////////////////////////////
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

    //! @brief ... get nonlocal radius
    //! @return ... nonlocal radius
    double GetNonlocalRadius() const;

    //! @brief ... set nonlocal radius
    //! @param rRadius...  nonlocal radius
    void SetNonlocalRadius(double rRadius);
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
    void Info(unsigned short rVerboseLevel) const;

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

    //! @brief ... nonlocal radius
    double mRadius;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check nonlocal radius is positive
    //! @param rRadius ... nonlocal radius
    void CheckNonlocalRadius(double rRadius) const;

};
}
#endif /* CONSTITUTIVENONLOCALDAMAGEPLASTICITY_H_ */
