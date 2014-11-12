/*
 * GradientDamageEngineeringStress.h
 *
 *  Created on: 10 Nov 2014
 *      Author: ttitsche
 */

#ifndef CONSTITUTIVEGRADIENTDAMAGEENGINEERINGSTRESS_H_
#define CONSTITUTIVEGRADIENTDAMAGEENGINEERINGSTRESS_H_

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class ConstitutiveStaticDataGradientDamagePlasticity1D;
class Logger;
class ConstitutiveTangentBase;
//! @author Joerg F. Unger
//! @date Apr 26, 2010
//! @brief ...
class GradientDamageEngineeringStress : public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    GradientDamageEngineeringStress();

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(
            ElementBase* rElement, int rIp,
            const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
            std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(
            ElementBase* rElement, int rIp,
            const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
            std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(
            ElementBase* rElement, int rIp,
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
    //! @brief ... get density
    //! @return ... density
    double GetDensity() const;

    //! @brief ... set density
    //! @param rRho ... density
    void SetDensity(double rRho);

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
    void SetNonlocalRadius(double rNonlocalRadius);

    //! @brief ... get thermal expansion coefficient
    //! @return ... thermal expansion coefficient
    double GetThermalExpansionCoefficient() const override;

    //! @brief ... set thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void SetThermalExpansionCoefficient(double rNu) override;

    //! @brief ... get damage law
    //! @return ... damage law
    NuTo::FullVector<double, Eigen::Dynamic> GetDamageLaw() const override;

    //! @brief ... set damage law
    //! @param rDamageLaw ... damage law
    void SetDamageLaw(const NuTo::FullVector<double, Eigen::Dynamic> rDamageLaw) override;

///////////////////////////////////////////////////////////////////////////

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

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const;

protected:
    //! @brief ... density
    double mRho;

    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... nonlocal radius
    double mNonlocalRadius;

    //! @brief ... thermal expansion coefficient \f$ \alpha \f$
    double mThermalExpansionCoefficient;

    //! @brief ... damage law type
    int mDamageLawType;

    //! @brief ... damage law parameters
    FullVector<double, Eigen::Dynamic> mDamageLawParameters;

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

    //! @brief ... check thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void CheckThermalExpansionCoefficient(double rAlpha) const;

    //! @brief ... check damage law
    //! @param rDamageLaw ... damage law
    void CheckDamageLaw(const NuTo::FullVector<double, Eigen::Dynamic>& rDamageLaw) const;

    //! @brief ... calculate the isotropic damage variable from the nonlocal eq strain history variable
    //! kappa and the damage law parameters (class members)
    //! \f$ \omega = 1 - \frac{\varepsilon_0}{\kappa} \exp \left(\frac{\varepsilon_0 - \kappa}{\varepsilon_f} \right) \f$
    //! @param rKappa ... history variable
    //! @return isotropic damage variable
    double CalculateDamage(double rKappa) const;

    //! @brief ... calculate the D_Omega_D_Kappa from the nonlocal eq strain history variable
    //! kappa and the damage law parameters (class members)
    //! \f$ \frac{\partial \omega}{\partial \kappa} = \frac{\varepsilon_0}{\kappa} \left(\frac{1}{\kappa} + \frac{1}{\varepsilon_f} \right) \exp \left(\frac{\varepsilon_0 - \kappa}{\varepsilon_f} \right) \f$
    //! @param rKappa ... history variable
    //! @return ... isotropic damage variable
    double CalculateDerivativeDamage(double rKappa) const;

};
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::GradientDamagePlasticityEngineeringStress)
#endif // ENABLE_SERIALIZATION

#endif /* CONSTITUTIVEGRADIENTDAMAGEENGINEERINGSTRESS_H_ */
