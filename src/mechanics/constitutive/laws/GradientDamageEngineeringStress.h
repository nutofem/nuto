#pragma once



#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"


namespace NuTo
{
namespace Constitutive
{
    enum class eDamageLawType;

}// namespace Constitutive

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
    typedef double StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<double>;

    GradientDamageEngineeringStress();

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLaw<GradientDamageEngineeringStress>>(*this, 0.0);
    }

    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
            const InterpolationType& rInterpolationType) const override;

    //! @brief Evaluate the constitutive relation in 2D
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param rStaticData Pointer to the history data.
    template<int TDim>
    eError Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                    const ConstitutiveOutputMap& rConstitutiveOutput,
                    Data& rStaticData);

    //! @brief Calculates the current static data based on the given CALCULATE_STATIC_DATA input.
    //! @param rStaticData History data.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @return Kappa value calculated from history data.
    double GetCurrentStaticData(Data& rStaticData, const ConstitutiveInputMap& rConstitutiveInput) const;

    //! @brief Calculates the error of the extrapolation.
    //! @param rStaticData History data.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @return Error of the extrapolation.
    double CalculateStaticDataExtrapolationError(Data& rStaticData, const ConstitutiveInputMap& rConstitutiveInput) const;

    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow,
                                                Node::eDof rDofCol,
                                                int rTimeDerivative) const override;

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

    //! @brief Evaluate the constitutive relation without handling the static data
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param rKappa new static data
    //! @param rKappaTangent new dKappa_dNonlocalEqStrain
    template<int TDim>
    eError EvaluateWithKappa(const ConstitutiveInputMap& rConstitutiveInput,
                             const ConstitutiveOutputMap& rConstitutiveOutput,
                             StaticDataType rKappa, double rKappaTangent);

    //! @brief Evaluate the static data part of the law
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param rStaticData Pointer to the history data.
    //! @return new static data kappa
    template<int TDim>
    double EvaluateStaticData(const ConstitutiveInputMap& rConstitutiveInput,
                              const ConstitutiveOutputMap& rConstitutiveOutput,
                              Data& rStaticData);


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
