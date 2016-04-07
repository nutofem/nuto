#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"


namespace NuTo
{
class ConstitutiveStaticDataGradientDamage;
template <int TDim> class EngineeringStrain;
template <int TRows, int TCols> class ConstitutiveMatrix;
template <int TRows> class ConstitutiveVector;
class ConstitutiveScalar;
class ConstitutiveStaticDataBase;


class Logger;
//! @author Thomas Titscher
//! @date Mar 14, 2015 // let's celebrate PI day!
//! @brief ...
class GradientDamageEngineeringStress : public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    GradientDamageEngineeringStress();

    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const override;

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override;

    //! @brief calculates the current static data based on the given CALCULATE_STATIC_DATA input
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @return ... current static data
    NuTo::ConstitutiveStaticDataGradientDamage GetCurrentStaticData(ElementBase& rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput) const;

    //! @brief calculates the error of the extrapolation
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @return ... error of the extrapolation
    double CalculateStaticDataExtrapolationError(ElementBase& rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput) const;


    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData1D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData2D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData3D(const ElementBase* rElement) const override;

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

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override {return false;}

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


protected:
    //! @brief ... density
    double mRho;

    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... nonlocal radius
    double mNonlocalRadius;

    //! @brief ... nonlocal radius parameter for non constant c-paramter
    double mNonlocalRadiusParameter;

    //! @brief ... thermal expansion coefficient \f$ \alpha \f$
    double mThermalExpansionCoefficient;

    //! @brief ... tensile strength
    double mTensileStrength;

    //! @brief ... uniaxial compressive strength
    double mCompressiveStrength;

    //! @brief ... fracture energy
    double mFractureEnergy;

    //! @brief ... damage law type
    Constitutive::eDamageLawType mDamageLawType;

private:

    //! @brief calculates the tangent \f$ \frac{1}{\xi} - \frac{\varepsilon_{eq} - \bar{\varepsilon}_{eq}}{\xi^2}\frac{\partial \xi}{\partial \varepsilon_{eq}}   \f$
    double CalculateLocalEqStrainXiFactor(double rLocalEqStrain, double rNonlocalEqStrain) const;

    //! @brief calculates \f$ \xi = c_0 + (c-c_0)\left(\frac{\varepsilon_{eq}}{e_{\xi}} \right) \f$
    double CalculateXi(double rLocalEqStrain) const;

};
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::GradientDamageEngineeringStress)
#endif // ENABLE_SERIALIZATION
