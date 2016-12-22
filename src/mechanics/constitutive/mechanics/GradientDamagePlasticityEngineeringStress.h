// $Id: GradientDamagePlasticityEngineeringStress.h 612 2012-08-13 07:31:23Z unger3 $
#ifndef CONSTITUTIVEGRADIENTDAMAGEPLASTICITYENGINEERINGSTRESS_H_
#define CONSTITUTIVEGRADIENTDAMAGEPLASTICITYENGINEERINGSTRESS_H_

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/elements/ElementEnum.h"

#include "mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class ConstitutiveStaticDataGradientDamagePlasticity1D;
class Logger;
class ConstitutiveTangentBase;
//! @author Joerg F. Unger
//! @date Apr 26, 2010
//! @brief ...
class GradientDamagePlasticityEngineeringStress : public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    GradientDamagePlasticityEngineeringStress();

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

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters()const override;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

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
    //! @param rStress              ... new stress
    //! @param rPlasticStrain       ... new plastic strain after return mapping
    //! @param rYieldConditionFlag  ... yield condition flag, true for active, false for inactive, 0 is Drucker-Prager, 1 is Rankine
    //! @param rDeltaKappa          ... delta equivalent plastic strain for Drucker Prager (0) and Rankine yield surface(1)
    //! @param rdSigmadEpsilon      ... new derivative of current stress with respect to the total strain
    //! @param rdKappadEpsilon1     ... new derivative of equivalent plastic strain (Drucker-Prager 0 Rankine 1) with respect to the total strain
    NuTo::Error::eError ReturnMapping1D(
            const EngineeringStrain1D& rStrain,
            const EngineeringStrain1D& rPrevPlasticStrain,
            const EngineeringStrain1D& rPrevTotalStrain,
            EngineeringStress1D& rStress,
            EngineeringStrain1D& rPlasticStrain,
            NuTo::FullMatrix<int,2,1>& rYieldConditionFlag,
            NuTo::FullMatrix<double,2,1>& rDeltaKappa,
            NuTo::FullMatrix<double,1,1>* rdSigma1dEpsilon1,
            NuTo::FullMatrix<double,2,1>* rdKappadEpsilon1,
            NuTo::Logger& rLogger)const;

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

    //! @brief ... performs the return mapping procedure for the plasticity model
    //! @param rStrain              ... current total strain
    //! @param rPrevPlasticStrain   ... previous plastic strain (history variable)
    //! @param rPrevTotalStrain     ... previous total strain (history variable)
    //! @param rPrevEqPlasticStrain ... previous equiavalente plastic strain (history variable)
    //! @param rEpsilonP            ... new plastic strain after return mapping
    //! @param rEqPlasticStrain     ... new equivalente olastic strain after return mapping
    //! @param rdSigmadEpsilon      ... new derivative of current stress with respect to the total strain
    //! @param rdEpsilonPdEpsilon   ... new derivative of current plastic strain with respect to the total strain
    NuTo::Error::eError ReturnMapping3D(
            const EngineeringStrain3D& rStrain,
            const double rPrevPlasticStrain[6],
            const EngineeringStrain3D& rPrevTotalStrain,
            Eigen::Matrix<double,6,1>& rStress,
            Eigen::Matrix<double,6,1>& rEpsilonP,
            double& rDeltaEqPlasticStrain,
            Eigen::Matrix<double,6,6>* rdSigmadEpsilon,
            Eigen::Matrix<double,6,6>* rdEpsilonPdEpsilon,
            NuTo::Logger& rLogger)const;


    //! @brief calculates the rankine yield surface and the derivatives with respect to the stress
    //! @param rStress current stress
    //! @param rFct tensile strength
    //! @param rdF_dSigma return value (first derivative)
    //! @param rd2F_d2Sigma return value (second derivative)
    //! @return yield function
    double YieldSurfaceRoundedRankine1D(Eigen::Matrix<double,1,1>& rStress, double rFct,
            Eigen::Matrix<double,1,1>* rdF_dSigma1,Eigen::Matrix<double,5,1>* rdF_dSigma2,
            Eigen::Matrix<double,1,1>* rd2F_d2Sigma1, Eigen::Matrix<double,5,1>* rd2F_dSigma2dSigma1)const;

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

    //! @brief calculates the Drucker Prager yield surface and the derivatives with respect to the stress
    //! @param rStress current stress
    //! @param rBETA parameter of the Drucker Prager yield surface
    //! @param rHP parameter of the Drucker Prager yield surface
    //! @param rdF_dSigma return value (first derivative)
    //! @param rd2F_d2Sigma return value (second derivative)
    //! @param rErrorDerivatives true, if derivative can't be calculated (on the hydrostatic axis)
    //! @return yield function
    double YieldSurfaceDruckerPrager1D(Eigen::Matrix<double,1,1>& rStress, double rBeta, double rHP,
            Eigen::Matrix<double,1,1>* rdF_dSigma1,Eigen::Matrix<double,5,1>* rdF_dSigma2,
            Eigen::Matrix<double,1,1>* rd2F_d2Sigma1, Eigen::Matrix<double,5,1>* rd2F_dSigma2dSigma1,
            bool &rErrorDerivatives
            )const;

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

    //! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
    //! @param rStress current stress
    //! @param rFct return tensile strength
    //! @param rdF_dSigma return value (first derivative)
    //! @param rd2F_d2Sigma return value (second derivative)
    //! @return yield function
    double YieldSurfaceRoundedRankine3D(
    		const Eigen::Matrix<double,6,1>& rStress,
    		double rFct,
    		Eigen::Matrix<double,6,1>* rdF_dSigma,
    		Eigen::Matrix<double,6,6>* rd2F_d2Sigma)const;

    //! @brief calculates the Drucker Prager yield surface and the derivatives with respect to the stress
    //! @param rStress current stress
    //! @param rBETA parameter of the Drucker Prager yield surface
    //! @param rHP parameter of the Drucker Prager yield surface
    //! @param rdF_dSigma return value (first derivative)
    //! @param rd2F_d2Sigma return value (second derivative)
    //! @param rErrorDerivatives true, if derivative can't be calculated (on the hydrostatic axis)
    //! @return yield function
    double YieldSurfaceDruckerPrager3D(Eigen::Matrix<double,6,1>& rStress, double rBeta, double rHP,
    		Eigen::Matrix<double,6,1>* rdF_dSigma,Eigen::Matrix<double,6,6>* rd2F_d2Sigma, bool &rErrorDerivatives
            )const;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override;

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

    //parameter related to the fracture energy
    FullMatrix<double,2,1> mEpsilonF;

    //! @brief ... yield mode (Drucker Prager, Rankine and combinations)
    Constitutive::eNonlocalDamageYieldSurface mYieldSurface;

    //! @brief ... overnonlocal parameter
    double mM;

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

    //! @brief ... calculate the isotropic damage variable from the nonlocal equivalente plastic strain
    //! @param rkappa scaled nonlocal equivalente plastic strain
    //! @param rKappaD material parameter related to the fracture energy and the nonlocal radius
    //! @return isotropic damage variable
    double CalculateDamage(double rNonlocalEpsilonPEq, double rKappaUnscaled)const;

    //! @brief ... calculate the isotropic damage variable from the local equivalente plastic strain and its derivative
    //! @param rkappa scaled local equivalente plastic strain
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
BOOST_CLASS_EXPORT_KEY(NuTo::GradientDamagePlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

#endif /* CONSTITUTIVEGRADIENTDAMAGEPLASTICITY_H_ */
