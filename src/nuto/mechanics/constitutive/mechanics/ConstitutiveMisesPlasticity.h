// $Id: LinearElastic.h 112 2009-11-17 16:43:15Z unger3 $

#ifndef ConstitutiveMisesPlasticity_H
#define ConstitutiveMisesPlasticity_H

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutivePiolaKirchhoffIIGreenLagrange.h"


namespace NuTo
{
class ConstitutiveStaticDataMisesPlasticity3D;

//! @brief ... mises plasticity with isotropic and kinematic hardening
//! @author Jörg F. Unger, ISM
//! @date December 2009
class ConstitutiveMisesPlasticity: public ConstitutiveEngineeringStressStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ConstitutiveMisesPlasticity();

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

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    void GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient1D& rDeformationGradient, double& rDamage) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    void GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient2D& rDeformationGradient, double& rDamage) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rDamage ... damage variable
    void GetDamage(const ElementBase* rElement, int rIp,
                                      const DeformationGradient3D& rDeformationGradient, double& rDamage) const;


    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient,
    		ConstitutiveTangentBase* rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    //! @param rTangent ... tangent
    void GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
    		const DeformationGradient3D& rDeformationGradient,
    		ConstitutiveTangentBase* rTangent) const;

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

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateTmpStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient1D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateTmpStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
    		const DeformationGradient2D& rDeformationGradient) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient
    void UpdateTmpStaticData_EngineeringStress_EngineeringStrain( ElementBase* rElement, int rIp,
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

    //! @brief ... performs the return mapping procedure in 3D
    //! @param rElement ... structure
    //! @param rIp ... integration point
    //! @param rDeformationGradient ... deformation gradient point
    //! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
    //! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
    //! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
    void ReturnMapping3D(const ElementBase* rElement,int rIp,
    		const DeformationGradient3D& rDeformationGradient,
    		EngineeringStress3D* rNewStress,
    		ConstitutiveTangentLocal6x6* rNewTangent,
    		ConstitutiveStaticDataMisesPlasticity3D* rNewStaticData)const;

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

    //! @brief ... get initial yield strength
    //! @return ... yield strength
    double GetInitialYieldStrength() const;

    //! @brief ... set initial yield strength
    //! @param rSigma ...  yield strength
    void SetInitialYieldStrength(double rSigma);

    //! @brief ... get yield strength for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding yield strength
    NuTo::FullMatrix<double> GetYieldStrength() const;

    //! @brief ... add yield strength
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  yield strength
    void AddYieldStrength(double rEpsilon, double rSigma);

    //! @brief ... get initial hardening modulus
    //! @return ... hardening modulus
    double GetInitialHardeningModulus() const;

    //! @brief ... set initial hardening modulus
    //! @param rH ...  hardening modulus
    void SetInitialHardeningModulus(double rH);

    //! @brief ... get hardening modulus for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding hardening modulus
    NuTo::FullMatrix<double> GetHardeningModulus() const;

    //! @brief ... add hardening modulus
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  hardening modulus
    void AddHardeningModulus(double rEpsilon, double rSigma);

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

    //! @brief ... equivalent strain with corresponding yield strength
    std::vector<std::pair<double, double> > mSigma;

    //! @brief ... equivalent strain with hardening modulus
    std::vector<std::pair<double, double> > mH;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check yield strength is positive
    //! @param rSigma ... yield strength
    void CheckYieldStrength(std::vector<std::pair<double, double> > rSigma) const;

    //! @brief ... check hardening modulus
    //! @param rH ... hardening modulus
    void CheckHardeningModulus(std::vector<std::pair<double, double> > rH) const;

    //! @brief ... calculates for a given equivalent plastic strain the radius of the yield surface
    //! @param ... rEpsilonPEq equivalent plastic strain
    //! @param ... rDSigmaDEpsilonP derivative of the yield strength with respect to the plastic strains (return value)
    //! @return ... yield strength (radius of the yield surface)
    double GetYieldStrength(double rEpsilonPEq, double& rDSigmaDEpsilonP)const;

    //! @brief ... calculates for a given equivalent plastic strain the hardening modulus
    //! @param ... rEpsilonPEq equivalent plastic strain
    //! @param ... rDSigmaDEpsilonP derivative of the hardening modulus with respect to the plastic strains (return value)
    //! @return ... hardening modulus
    double GetHardeningModulus(double rEpsilonPEq, double& rDHDEpsilonP)const;
};
}  //namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveMisesPlasticity)

//this is due to the diamond structure (virtual public derivative)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::ConstitutiveEngineeringStressStrain, NuTo::ConstitutiveMisesPlasticity>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION
#endif // ConstitutiveMisesPlasticity_H_
