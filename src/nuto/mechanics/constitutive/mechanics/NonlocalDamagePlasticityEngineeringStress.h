// $Id$ 
#ifndef CONSTITUTIVENONLOCALDAMAGEPLASTICITYENGINEERINGSTRESS_H_
#define CONSTITUTIVENONLOCALDAMAGEPLASTICITYENGINEERINGSTRESS_H_

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain;
class Logger;
class ConstitutiveTangentBase;
//! @author Joerg F. Unger
//! @date Apr 26, 2010
//! @brief ...
class NonlocalDamagePlasticityEngineeringStress : public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    NonlocalDamagePlasticityEngineeringStress();

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
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const override;


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
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief ... performs the return mapping procedure for the plasticity model
    //! @param rStrain              ... current total strain
    //! @param rPrevPlasticStrain   ... previous plastic strain (history variable)
    //! @param rPrevTotalStrain     ... previous total strain (history variable)
    //! @param rPrevEqPlasticStrain ... previous equiavalente plastic strain (history variable)
    //! @param rEpsilonP            ... new plastic strain after return mapping
    //! @param rEqPlasticStrain     ... new equivalente olastic strain after return mapping
    //! @param rdEpsilonPdEpsilon   ... new derivative of current plastic strain with respect to the total strain
    NuTo::Error::eError ReturnMapping2D(
            const EngineeringStrain2D& rStrain,
            const double rPrevPlasticStrain[4],
            const EngineeringStrain2D& rPrevTotalStrain,
            Eigen::Matrix<double,4,1>& rStress,
            Eigen::Matrix<double,4,1>& rEpsilonP,
                double& rDeltaEqPlasticStrain,
            Eigen::Matrix<double,4,4>& rdEpsilonPdEpsilon,
            Logger& rLogger)const;

    //! @brief calculates the rounded rankine yield surface
    //! @param rStress current stress
    //! @param rSigma_1 first principal stress
    //! @param rSigma_2 second principal stress
    //! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
    //! @return yield condition
    double YieldSurfaceRankine2DRounded(Eigen::Matrix<double,4,1>& rStress, double rFct)const;

    //! @brief calculates the first and second derivative of the rounded Rankine yield surface with respect to the stress
    //! @param dF_dsigma return value (first derivative)
    //! @param d2F_d2sigma return value (second derivative)
    //! @param stress vector
    //! @param sigma_1 first principal stress in the plane
    //! @param sigma_2 second principal stress in the plane
    //! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm value_sqrt$
    void YieldSurfaceRankine2DRoundedDerivatives(Eigen::Matrix<double,4,1>& rdF_dSigma,Eigen::Matrix<double,4,4>* rd2F_d2Sigma,
    		Eigen::Matrix<double,4,1>& rStress)const;

    //! @brief calculates the drucker prager yield surface
    //! @param rStress current stress
    //! @param rBeta parameter of the Drucker-Prager yield surface
    //! @param rHP parameter of the Drucker-Prager yield surface
    //! @return yield condition
    double YieldSurfaceDruckerPrager2D(Eigen::Matrix<double,4,1>& rStress, double rBeta, double rHP)const;

    //! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
    //! @param dF_dsigma return value (first derivative)
    //! @param d2F_d2sigma return value (second derivative)
    //! @param elasticStress current stress
     //! @param BETA parameter of the Drucker Prager yield surface
    //! @return false, if the stress is on the hydrostatic axis otherwise true
    bool YieldSurfaceDruckerPrager2DDerivatives(Eigen::Matrix<double,4,1>& dF_dsigma,Eigen::Matrix<double,4,4>* d2F_d2sigma,
    		Eigen::Matrix<double,4,1>& elasticStress, double BETA)const;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const;

    //! @brief ... returns true, if a material model has is nonlocal (stiffness is of dynamic size, nonlocal averaging)
    //! @return ... see brief explanation
    bool IsNonlocalModel()const;

protected:
    //! @brief ... density
    double mRho;

    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... nonlocal radius
    double mNonlocalRadius;

    //! @brief ... uniaxial tensile strength
    double mTensileStrength;

    //! @brief ... uniaxial compressive strength
    double mCompressiveStrength;

    //! @brief ... biaxial compressive strength
    double mBiaxialCompressiveStrength;

    //! @brief ... fracture energy
    double mFractureEnergy;

    //! @brief ... yield mode (Drucker Prager, Rankine and combinations)
    Constitutive::eNonlocalDamageYieldSurface mYieldSurface;

    //! @brief ... either use a pure plasticity model (false) or add softening using the damage model (true)
    bool mDamage;

    //! @brief ... thermal expansion coefficient \f$ \alpha \f$
    double mThermalExpansionCoefficient;

    //! @brief ... check if density is non negativ
    //! @param rRho ... Young's modulus
    void CheckDensity(double rRho) const;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check if nonlocal radius is positive
    //! @param rRadius ... nonlocal radius
    void CheckNonlocalRadius(double rNonlocalRadius) const;

    //! @brief ... check if tensile strength is positive
    //! @param rTensileStrength ... nonlocal radius
    void CheckTensileStrength(double rTensileStrength) const;

    //! @brief ... check if compressive strength is positive
    //! @param rRadius ... compressive strength
    void CheckCompressiveStrength(double rCompressiveStrength) const;

    //! @brief ... check if biaxial compressive strength is positive
    //! @param rBiaxialCompressiveStrength ... biaxial compressive strength
    void CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const;

    //! @brief ... check if fracture energy is positive
    //! @param rFractureEnergy ... fracture energy
    void CheckFractureEnergy(double rFractureEnergy) const;

    //! @brief ... check thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void CheckThermalExpansionCoefficient(double rAlpha) const;

    //! @brief ... calculate the length of the element in plane coordinates (square root of area)
    //! @brief ... an interpolation is made based on the principal elastic stress in plane and thickness direction
    double CalculateEquivalentLength2D(const ElementBase* rElement, const Eigen::Matrix<double,4,1>& rStress) const;

    //! @brief ... calculates the derivative of the equivalent length of the element in plane coordinates (square root of area) with respect to the plastic strains
    //! @brief ... an interpolation is made based on the principal strains in plane and thickness direction
    double CalculateDerivativeEquivalentLength2D(const ElementBase* rElement, const Eigen::Matrix<double,4,1>& rStress,
    		                    Eigen::Matrix<double,4,1>& rdLdStress) const;

    //! @brief ... calculate the nonlocal equivalente plastic strain of an integration point
    //! @param rElement Element
    //! @param rIp integration point
    //! @param rUnloading true, if all nonlocal integration points are in an unloading state
    //! @return equivalente plastic strain
    double CalculateNonlocalEquivalentPlasticStrain(const ElementBase* rElement, int rIp, bool& rUnloading)const;

    //! @brief ... calculate the nonlocal plastic strain of an integration point
    //! @param rElement Element
    //! @param rIp integration point
    //! @param rNonlocalPlasticStrain nonlocal plastic strain (return value)
    //! @return void
    void CalculateNonlocalPlasticStrain(const ElementBase* rElement, int rIp, double rNonlocalPlasticStrain[4])const;

    //! @brief ... calculate the isotropic damage variable from the nonlocal equivalente plastic strain
    //! @param rkappa scaled nonlocal equivalente plastic strain
    //! @param rKappaD material parameter related to the fracture energy and the nonlocal radius
    //! @return isotropic damage variable
    double CalculateDamage(double rNonlocalEpsilonPEq, double rKappaUnscaled)const;

    //! @brief ... calculate the isotropic damage variable from the nonlocal equivalente plastic strain and its derivative
    //! @param rkappa scaled nonlocal equivalente plastic strain
    //! @param rKappaD material parameter related to the fracture energy and the nonlocal radius
    //! @param rdOmegadKappa derivative of omega w.r.t. Kappa
    //! @return isotropic damage variable
    double CalculateDerivativeDamage(double rkappa, double rKappaD, double& rdOmegadKappa)const;

    //! @brief ... calculate the unscaled kappa_d used for scaling the equivalent plastic strain (related to the fracture energy)
    //! @brief ... return kappa_d
    double CalculateKappaD()const;
};
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NonlocalDamagePlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

#endif /* CONSTITUTIVENONLOCALDAMAGEPLASTICITY_H_ */
