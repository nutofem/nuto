// $Id$

#ifndef LINEARELASTICENGINEERINGSTRESS_H_
#define LINEARELASTICENGINEERINGSTRESS_H_

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
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
//! @date July 2012
class LinearElasticEngineeringStress: public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    LinearElasticEngineeringStress();

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

    ///////////////////////////////////////////////////////////////////////////

    // calculate coefficients of the material matrix
    void CalculateCoefficients2DPlainStress(double& C11, double& C12, double& C33) const;
    void CalculateCoefficients3D(double& C11, double& C12, double& C44) const;

    // parameters /////////////////////////////////////////////////////////////
    //! @brief ... get density
    //! @return ... density
    virtual double GetDensity() const override;

    //! @brief ... set density
    //! @param rRho ... density
    virtual void SetDensity(double rRho) override;

    //! @brief ... get Young's modulus
    //! @return ... Young's modulus
    double GetYoungsModulus() const override;

    //! @brief ... set Young's modulus
    //! @param rE ... Young's modulus
    void SetYoungsModulus(double rE) override;

    //! @brief ... get Poisson's ratio
    //! @return ... Poisson's ratio
    double GetPoissonsRatio() const override;

    //! @brief ... set Poisson's ratio
    //! @param rNu ... Poisson's ratio
    void SetPoissonsRatio(double rNu) override;

    //! @brief ... get thermal expansion coefficient
    //! @return ... thermal expansion coefficient
    double GetThermalExpansionCoefficient() const override;

    //! @brief ... set thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void SetThermalExpansionCoefficient(double rNu) override;
    ///////////////////////////////////////////////////////////////////////////

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters()const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    //! @param rLogger stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override
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
BOOST_CLASS_EXPORT_KEY(NuTo::LinearElasticEngineeringStress)
#endif //ENABLE_SERIALIZATION

#endif // LINEARELASTICENGINEERINGSTRESS_H_
