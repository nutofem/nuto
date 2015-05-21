// $Id$

#ifndef CONSTITUTIVEBASE_H_
#define CONSTITUTIVEBASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <string>
#include <vector>
#include <map>

#include "nuto/base/ErrorEnum.h"
#include "nuto/math/FullMatrix_Def.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
// forward declarations
class ConstitutiveEngineeringStressStrain;
class ConstitutiveLatticeStressStrain;
class ConstitutiveStaticDataBase;
class ConstitutiveTangentLocal1x1;
class ConstitutiveTangentLocal2x2;
class ConstitutiveTangentLocal3x3;
class ConstitutiveTangentLocal6x6;
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class ElementBase;
class EngineeringStrain3D;
class EngineeringStress1D;
class EngineeringStress2D;
class EngineeringStress3D;
class Logger;
class SecondPiolaKirchhoffStress1D;
class SecondPiolaKirchhoffStress2D;
class SecondPiolaKirchhoffStress3D;
class StructureBase;

//! @brief ... base class for the constitutive relationship, e.g. material laws
//! @author Stefan Eckardt, ISM
//! @date November 2009
class ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    ConstitutiveBase();

    //! @brief ... constructor
    virtual ~ConstitutiveBase()
    {}

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput);

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput);

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput);

    // parameters /////////////////////////////////////////////////////////////

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    virtual bool GetVariableBool(Constitutive::eConstitutiveVariable rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    virtual void SetVariableBool(Constitutive::eConstitutiveVariable rIdentifier, bool rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    virtual double GetVariableDouble(Constitutive::eConstitutiveVariable rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    virtual void SetVariableDouble(Constitutive::eConstitutiveVariable rIdentifier, double rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    virtual NuTo::FullVector<double,Eigen::Dynamic> GetVariableFullVectorDouble(Constitutive::eConstitutiveVariable rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    virtual void SetVariableFullVectorDouble(Constitutive::eConstitutiveVariable rIdentifier, NuTo::FullVector<double,Eigen::Dynamic> rValue);


    //! @brief ... get density
    //! @return ... density
    virtual double GetDensity() const;

    //! @brief ... set density
    //! @param rRho ... density
    virtual void SetDensity(double rRho);

    //! @brief ... get Young's modulus
    //! @return ... Young's modulus
    virtual double GetYoungsModulus() const;

    //! @brief ... get factor to modify Youngs modulus (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorYoungsModulus(const ElementBase* rElement,int rIp) const;

    //! @brief ... set Young's modulus
    //! @param rE ... Young's modulus
    virtual void SetYoungsModulus(double rE);

    //! @brief ... get Poisson's ratio
    //! @return ... Poisson's ratio
    virtual double GetPoissonsRatio() const;

    //! @brief ... set Poisson's ratio
    //! @param rNu ... Poisson's ratio
    virtual void SetPoissonsRatio(double rNu);

    //! @brief ... get thermal expansion coefficient
    //! @return ... thermal expansion coefficient
    virtual double GetThermalExpansionCoefficient() const;

    //! @brief ... set thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    virtual void SetThermalExpansionCoefficient(double rAlpha);

    //! @brief ... get factor to modify Poisson's ratio (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorPoissonsRatio(const ElementBase* rElement,int rIp) const;

    //! @brief ... get initial yield strength
    //! @return ... yield strength
    virtual double GetInitialYieldStrength() const;

    //! @brief ... set initial yield strength
    //! @param rSigma ...  yield strength
    virtual void SetInitialYieldStrength(double rSigma);

    //! @brief ... get yield strength for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding yield strength
    virtual NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetYieldStrength() const;

    //! @brief ... add yield strength
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  yield strength
    virtual void AddYieldStrength(double rEpsilon, double rSigma);

    //! @brief ... get factor to modify yield strength (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorYieldStrength(const ElementBase* rElement,int rIp) const;

    //! @brief ... get initial hardening modulus
    //! @return ... hardening modulus
    virtual double GetInitialHardeningModulus() const;

    //! @brief ... set hardening modulus
    //! @param rH ...  hardening modulus
    virtual void SetInitialHardeningModulus(double rH);

    //! @brief ... get hardening value
    //! @return ... hardening value
    virtual double GetHardeningValue() const;

    //! @brief ... set hardening value
    //! @param rH ...  hardening value
    virtual void SetHardeningValue(double rH);

    //! @brief ... get hardening exponent
    //! @return ... hardening exponent
    virtual double GetHardeningExponent() const;

    //! @brief ... set hardening exponent
    //! @param rHexponent ...  hardening exponent
    virtual void SetHardeningExponent(double rHexponent);

    //! @brief ... get hardening modulus for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding hardening modulus
    virtual NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetHardeningModulus() const;

    //! @brief ... add hardening modulus
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  hardening modulus
    virtual void AddHardeningModulus(double rEpsilon, double rH);

    //! @brief ... get factor to modify hardening modulus (using random fields)
    //! @param rElement ...  element
    //! @param rIp ...  integration point
    double GetRanfieldFactorHardeningModulus(const ElementBase* rElement,int rIp) const;

    //! @brief ... get nonlocal radius
    //! @return ... nonlocal radius
    virtual double GetNonlocalRadius() const;

    //! @brief ... set nonlocal radius
    //! @param rRadius...  nonlocal radius
    virtual void SetNonlocalRadius(double rRadius);

    //! @brief ... get nonlocal radius parameter
    //! @return ... nonlocal radius parameter
    virtual double GetNonlocalRadiusParameter() const;

    //! @brief ... set nonlocal radius parameter
    //! @param rRadius...  nonlocal radius parameter
    virtual void SetNonlocalRadiusParameter(double rRadiusParameter);

    //! @brief ... get tensile strength
    //! @return ... tensile strength
    virtual double GetTensileStrength() const;

    //! @brief ... set tensile strength
    //! @param rTensileStrength...  tensile strength
    virtual void SetTensileStrength(double rTensileStrength);

    //! @brief ... get shear strength
    //! @return ... shear strength
    virtual double GetShearStrength() const;

    //! @brief ... set shear strength
    //! @param rShearStrength...  shear strength
    virtual void SetShearStrength(double rShearStrength);

    //! @brief ... get compressive strength
    //! @return ... compressive strength
    virtual double GetCompressiveStrength() const;

    //! @brief ... set compressive strength
    //! @param rCompressiveStrength...  compressive strength
    virtual void SetCompressiveStrength(double rCompressiveStrength);

    //! @brief ... get biaxial compressive strength
    //! @return ... biaxial compressive strength
    virtual double GetBiaxialCompressiveStrength() const;

    //! @brief ... set biaxial compressive strength
    //! @param rBiaxialCompressiveStrength...  biaxial compressive strength
    virtual void SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength);

    //! @brief ... get fracture energy
    //! @return ... fracture energy
    virtual double GetFractureEnergy() const;

    //! @brief ... set fracture energy
    //! @param rFractureEnergy... fracture energy
    virtual void SetFractureEnergy(double rFractureEnergy);

    //! @brief ... get friction coefficient
    //! @return ... friction coefficient
    virtual double GetFrictionCoefficient() const;

    //! @brief ... set friction coefficient
    //! @param rFrictionCoefficient... friction coefficient
    virtual void SetFrictionCoefficient(double rFrictionCoefficient);

    //! @brief ... get HeatCapacity
    //! @return ... HeatCapacity
    virtual double GetHeatCapacity() const;

    //! @brief ... set HeatCapacity
    //! @param rHeatCapacity ... HeatCapacity
    virtual void SetHeatCapacity(double rHeatCapacity);

    //! @brief ... get thermal conductivity
    //! @return ... thermal conductivity
    virtual double GetThermalConductivity() const;

    //! @brief ... set thermal conductivity
    //! @param ... thermal conductivity
    virtual void SetThermalConductivity(double rThermalConductivity);

    //! @brief ... get viscosity
    //! @return ... viscosity
    virtual double GetViscosity() const;

    //! @brief ... set viscosity
    //! @param ... viscosity
    virtual void SetViscosity(double rViscosity);

    //! @brief ... get viscosity exponent
    //! @return ... viscosity exponent
    virtual double GetViscosityExponent() const;

    //! @brief ... set viscosity exponent
    //! @param ... viscosity exponent
    virtual void SetViscosityExponent(double rViscosityExponent);

    //! @brief ... get damage distribution (determines the portion of damage via viscoplasticity and plasticity)
    //! @return ... damage distribution
    virtual double GetDamageDistribution() const;

    //! @brief ... set damage distribution (determines the portion of damage via viscoplasticity and plasticity)
    //! @param ... damage distribution
    virtual void SetDamageDistribution(double rDamageDistribution);

    //! @brief ... get damage law
    //! @return ... damage law
    virtual NuTo::FullVector<double, Eigen::Dynamic> GetDamageLaw() const;

    //! @brief ... set damage law
    //! @param rDamageLawParameters ... damage law
    virtual void SetDamageLaw(const NuTo::FullVector<double, Eigen::Dynamic> rDamageLaw);

    //! @brief ... get viscoplastic yield surface offset with respect to the plastic yield surface
    //! @return ... viscoplastic yield surface offset
    virtual double GetViscoplasticYieldSurfaceOffset() const;

    //! @brief ... set viscoplastic yield surface offset with respect to the plastic yield surface
    //! @param ... viscoplastic yield surface offset
    virtual void SetViscoplasticYieldSurfaceOffset(double rViscoplasticYieldSurfaceOffset);

    //! @brief ... get fatigue extrapolation flag
    //! @param ... fatigue extrapolation flag
    virtual bool GetFatigueExtrapolation() const;

    //! @brief ... set fatigue extrapolation flag
    //! @param ... fatigue extrapolation flag
    virtual void SetFatigueExtrapolation(bool rFatigueExtrapolation);

    //! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
    //! @param rRelativeHumidity ... relative humidity
    //! @return ... equilibrium water volume fraction
    virtual double GetEquilibriumWaterVolumeFraction(double rRelativeHumidity, NuTo::FullVector<double,Eigen::Dynamic> rCoeffs) const;


    ///////////////////////////////////////////////////////////////////////////

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const = 0;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const = 0;

    //! @brief ... returns whether the parameters of the constitutive relationship are valid or not
    //! @return ...  <B>true</B> if all parameters of the constitutive relationship are valid and <B>false</B> otherwise
    inline bool AreParametersValid() const
    {
        return this->mParametersValid;
    }

    //! @brief ... check if all parameters are valid and modify parameter validity flag
    //! @sa mParametersValid
    void SetParametersValid();

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel, Logger& rLogger) const;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const=0;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual const ConstitutiveEngineeringStressStrain* AsConstitutiveEngineeringStressStrain()const;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual ConstitutiveEngineeringStressStrain* AsConstitutiveEngineeringStressStrain();

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual const ConstitutiveLatticeStressStrain* AsConstitutiveLatticeStressStrain()const;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual ConstitutiveLatticeStressStrain* AsConstitutiveLatticeStressStrain();

    //! @brief ... allocate the correct static data
    //! @return ... see brief explanation
    virtual ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement)const;

    //! @brief ... allocate the correct static data
    //! @return ... see brief explanation
    virtual ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement)const;

    //! @brief ... allocate the correct static data
    //! @return ... see brief explanation
    virtual ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void CheckParameters() const = 0;

    //! @brief ... flag which is <B>true</B> if all parameters of the constitutive relationship are valid and <B>false</B> otherwise
    bool mParametersValid;


};

}

#endif // CONSTITUTIVEBASE_H_ 
