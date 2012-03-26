// $Id:$

#ifndef CONSTITUTIVE_LATTICE_CONCRETE_H
#define CONSTITUTIVE_LATTICE_CONCRETE_H

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveLatticeStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutivePiolaKirchhoffIIGreenLagrange.h"

namespace NuTo
{

//! @brief ... lattice model to simulate concrete
//! @author Jörg F. Unger, NU
//! @date January 2012
class ConstitutiveLatticeConcrete: public ConstitutiveLatticeStressStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ConstitutiveLatticeConcrete();

    //  Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rLatticeStrain ... Lattice strain
    Error::eError GetLatticePlasticStrain(const ElementBase* rElement, int rIp,
                                      const LatticeStrain2D& rLatticeStrain, LatticeStrain3D& rLatticePlasticStrain) const;

    //  Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice plastic strain from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rLatticeStrain ... Lattice strain
    Error::eError GetLatticePlasticStrain(const ElementBase* rElement, int rIp,
                                      const LatticeStrain3D& rLatticeStrain, LatticeStrain3D& rLatticePlasticStrain) const;

    // Lattice stress - Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rLatticeStress ... Lattice stress
    Error::eError GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
    		      const LatticeStrain2D& rLatticeStrain, LatticeStress2D& rLatticeStress) const;

    // Lattice stress - Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rLatticeStress ... Lattice stress
    Error::eError GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
    		      const LatticeStrain2D& rLatticeStrain, LatticeStress3D& rLatticeStress) const;

    // Lattice stress - Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the deformation gradient)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rCauchyStress ... Cauchy stress
    Error::eError GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
    		      const LatticeStrain3D& rLatticeStrain, LatticeStress3D& rLatticeStress) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rDamage ... damage variable
    Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const LatticeStrain2D& rLatticeStrain, double& rDamage) const;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from deformation gradient in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rDamage ... damage variable
    Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const LatticeStrain3D& rLatticeStrain, double& rDamage) const;

    //! @brief ... calculate the tangent (derivative of the Lattice stresses with respect to the Lattice strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
    		const LatticeStrain2D& rLatticeStrain,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... calculate the tangent (derivative of the Lattice stresses with respect to the Lattice strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    //! @param rTangent ... tangent
    Error::eError GetTangent_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
    		const LatticeStrain3D& rLatticeStrain,
            ConstitutiveTangentBase* rTangent) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError UpdateStaticData_LatticeStress_LatticeStrain( ElementBase* rElement, int rIp,
    		const LatticeStrain2D& rLatticeStrain) const;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError UpdateStaticData_LatticeStress_LatticeStrain( ElementBase* rElement, int rIp,
    		const LatticeStrain3D& rLatticeStrain) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError UpdateTmpStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain) const;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError UpdateTmpStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataLatticeStress_LatticeStrain2D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataLatticeStress_LatticeStrain3D(const ElementBase* rElement) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError GetInternalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain, double& rEnergy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError GetInternalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain, double& rEnergy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... deformation gradient
    Error::eError GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain, double& rEnergy) const;
    ///////////////////////////////////////////////////////////////////////////

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

    //! @brief ... get tensile strength
    //! @return ... tensile strength
    double GetTensileStrength() const;

    //! @brief ... set tensile strength
    //! @param rTensileStrength...  tensile strength
    void SetTensileStrength(double rTensileStrength);

    //! @brief ... get shear strength
    //! @return ... shear strength
    double GetShearStrength() const;

    //! @brief ... set shear strength
    //! @param rShearStrength...  shear strength
    void SetShearStrength(double rShearStrength);

    //! @brief ... get fracture energy
    //! @return ... fracture energy
    double GetFractureEnergy() const;

    //! @brief ... set fracture energy
    //! @param rFractureEnergy... fracture energy
    void SetFractureEnergy(double rFractureEnergy);

    //! @brief ... get friction coefficient
    //! @return ... friction coefficient
    double GetFrictionCoefficient() const;

    //! @brief ... set friction coefficient
    //! @param rFrictionCoefficient... friction coefficient
    void SetFrictionCoefficient(double rFrictionCoefficient);
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

    //! @brief ... shear factor \f$ \alpha \f$
    double mAlpha;

    //! @brief ... shear strength \f$ \sigma_t \f$
    double mSigmaT;

    //! @brief ... normal strength \f$ \sigma_n \f$
    double mSigmaN;

    //! @brief ... fracture energy \f$ G \f$
    double mG;

    //! @brief ...  \f$ n_t \f$
    double mNt;

    //! @brief ... friction coefficient \f$ \mu \f$
    double mMu;

    //! @brief ... check if density is positive
    //! @param rRho ... density
    void CheckDensity(double rRho) const;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check if tensile strength is positive
    //! @param rTensileStrength ... tensile strength
    void CheckTensileStrength(double rTensileStrength) const;

    //! @brief ... check if shear strength is positive
    //! @param rShearStrength ... shear strength
    void CheckShearStrength(double rShearStrength) const;

    //! @brief ... check if fracture energy is positive
    //! @param rFractureEnergy ... fracture energy
    void CheckFractureEnergy(double rFractureEnergy) const;
    //! @brief ... check friction coefficient is positive
    //! @param rFrictionCoefficient ... friction coefficient
    void CheckFrictionCoefficient(double rFrictionCoefficient) const;

    //! @brief ... calculate the elastic parameters for the plane
    void CalculateElasticParameters(double& rEn,double& rEt)const;

    //! @brief ... calculate the stress and stiffness in tension (positive normal strain)
    //! param rElement element
    //! param rIp integration point
    //! param rLatticeStrain lattice strain
    //! param rLatticeStress lattice stress
    void CutOffTension2D(const ElementBase* rElement, int rIp,
    		  const LatticeStrain2D& rLatticeStrain, LatticeStress2D *rLatticeStress, ConstitutiveTangentLocal2x2 *rTangent)const;

    //! @brief ... calculate the stress and  in compression (negative normal strain)
    //! param rElement element
    //! param rIp integration point
    //! param rLatticeStrain lattice strain
    //! param rLatticeStress lattice stress
    void CutOffCompression2D(const ElementBase* rElement, int rIp,
  		  const LatticeStrain2D& rLatticeStrain, LatticeStress2D *rLatticeStress, ConstitutiveTangentLocal2x2 *rTangent)const;


 };

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveLatticeConcrete)

//this is due to the diamond structure (virtual public derivative)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::ConstitutiveLatticeStressStrain, NuTo::ConstitutiveLatticeConcrete>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVE_LATTICE_CONCRETE_H
