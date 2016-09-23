//============================================================================
// Name        : LocalDamageModel.cpp
// Author      : Philip Huschke
// Version     : 01 Sep 2016
// Copyright   :
// Description : Local damage model
//============================================================================

#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{
template <int TRows, int TCols> class ConstitutiveMatrix;
template <int TRows> class ConstitutiveVector;
class ConstitutiveScalar;
class ConstitutiveStaticDataBase;
class ConstitutiveStaticDataGradientDamage;
enum class eSectionType;
class LocalDamageModel: public ConstitutiveBase
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    LocalDamageModel();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


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
    NuTo::eError Evaluate1D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::eError Evaluate2D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override;

    double Evaluate2D(NuTo::eSectionType rSectionType, const ConstitutiveStaticDataGradientDamage& rCurrentStaticData, const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::eError Evaluate3D(ElementBase* rElement, int rIp,
            const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput) override;


    //! @brief calculates the current static data based on the given CALCULATE_STATIC_DATA input
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @return ... current static data
    ConstitutiveStaticDataGradientDamage GetCurrentStaticData(ElementBase& rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput) const;

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters() const override;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData1D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData2D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData3D(const ElementBase* rElement) const override;

    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow,
                                                Node::eDof rDofCol,
                                                int rTimeDerivative) const override;

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

private:
    //! @brief Calculates the damage
    double CalculateDamage(double rKappa) const;

    double CalculateDerivativeDamage(double rKappa) const;

    //! @brief Density \f$ \rho \f$
    double mDensity;

    //! @brief Young's modulus \f$ E \f$
    double mYoungsModulus;

    //! @brief Poisson's ratio \f$ \nu \f$
    double mPoissonsRatio;

    //! @brief Tensile strength \f$ \sigma_t \f$
    double mTensileStrength;

    //! @brief Compressive strength \f$ \sigma_c \f$
    double mCompressiveStrength;

    //! @brief Fracture energy \f$ G_c \f$
    double mFractureEnergy;

    //! @brief Damage law type
    Constitutive::eDamageLawType mDamageLawType;

    //! @brief Maximum damage
    static constexpr double mMaxDamage = 0.90;

    //! @brief Parameter for the exponential damage function
    static constexpr double mAlpha = 0.99;

};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LocalDamageModel)
#endif // ENABLE_SERIALIZATION
