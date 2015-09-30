// $Id: LinearElastic.h 112 2009-11-17 16:43:15Z unger3 $

#ifndef MisesPlasticityEngineeringStress_H
#define MisesPlasticityEngineeringStress_H

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
namespace NuTo
{
class ConstitutiveStaticDataMisesPlasticity3D;
class Logger;
template <int TNumRows,int TNumColumns>
class ConstitutiveTangentLocal;
//! @brief ... mises plasticity with isotropic and kinematic hardening
//! @author Jörg F. Unger, ISM
//! @date December 2009
class MisesPlasticityEngineeringStress: public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    MisesPlasticityEngineeringStress();

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const;

    //! @brief ... performs the return mapping procedure in 2D
    //! @param rElement ... structure
    //! @param rIp ... integration point
    //! @param rEngineeringStrain ... engineering strain
    //! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
    //! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
    //! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
    NuTo::Error::eError ReturnMapping2D(const ElementBase* rElement,int rIp,
            const EngineeringStrain2D& rEngineeringStrain,
            EngineeringStress2D* rNewStress,
            ConstitutiveTangentLocal<3,3>* rNewTangent,
            ConstitutiveStaticDataMisesPlasticity3D* rNewStaticData,
            Logger& rLogger)const;

    //! @brief ... performs the return mapping procedure in 3D
    //! @param rElement ... structure
    //! @param rIp ... integration point
    //! @param rEngineeringStrain ... engineering strain
    //! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
    //! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
    //! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
    NuTo::Error::eError ReturnMapping3D(const ElementBase* rElement,int rIp,
    		const EngineeringStrain3D& rEngineeringStrain,
    		EngineeringStress3D* rNewStress,
    		ConstitutiveTangentLocal<6,6>* rNewTangent,
    		ConstitutiveStaticDataMisesPlasticity3D* rNewStaticData,
    		Logger& rLogger)const;

    // calculate coefficients of the material matrix
    void CalculateCoefficients3D(double& C11, double& C12, double& C44) const;

    // parameters /////////////////////////////////////////////////////////////

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... get yield strength for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding yield strength
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetYieldStrength() const override;

    //! @brief ... add yield strength
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  yield strength
    void AddYieldStrength(double rEpsilon, double rSigma) override;

    //! @brief ... get hardening modulus for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding hardening modulus
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetHardeningModulus() const  override;

    //! @brief ... add hardening modulus
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  hardening modulus
    void AddHardeningModulus(double rEpsilon, double rSigma) override;

    //! @brief ... get dimension of the constitutive relationship
    //! @return ... dimension of the constitutive relationship (1, 2 or 3)
    int GetGlobalDimension() const;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const  override;

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

    //! @brief ... energy flag (store prev stress and strain to calculate energy using trapezoidal rule
    bool mEnergyFlag;

    //! @brief ... equivalent strain with corresponding yield strength
    std::vector<std::pair<double, double> > mSigma;

    //! @brief ... equivalent strain with hardening modulus
    std::vector<std::pair<double, double> > mH;

    //! @brief ... thermal expansion coefficient
    double mThermalExpansionCoefficient;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check yield strength is positive
    //! @param rSigma ... yield strength
    void CheckYieldStrength(std::vector<std::pair<double, double> > rSigma) const;

    //! @brief ... check thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void CheckThermalExpansionCoefficient(double rAlpha) const;

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
BOOST_CLASS_EXPORT_KEY(NuTo::MisesPlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION
#endif // ConstitutiveMisesPlasticity_H_
