//
// Created by Thomas Titscher on 10/24/16.
//
#pragma once
#include "nuto/mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
namespace NuTo
{

class GradientDamageFatigueEngineeringStress : public GradientDamageEngineeringStress
{
public:

    typedef Eigen::Vector2d StaticDataType;
    using Data = typename Constitutive::IPConstitutiveLaw<GradientDamageFatigueEngineeringStress>::Data;

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLaw<GradientDamageFatigueEngineeringStress>>(*this, Eigen::Vector2d(0.,0.));
    }

    //! @brief Evaluate the constitutive relation in 2D
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param rStaticData Pointer to the history data.
    template<int TDim>
    eError Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                    const ConstitutiveOutputMap& rConstitutiveOutput,
                    Data& rStaticData)
    {
        double kappa = GetCurrentStaticData(rStaticData, rConstitutiveInput);
        return GradientDamageEngineeringStress::EvaluateWithKappa<TDim>(rConstitutiveInput, rConstitutiveOutput, kappa);
    }

    //! @brief Calculates the current static data based on the given CALCULATE_STATIC_DATA input.
    //! @param rStaticData History data.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @return Kappa value calculated from history data.
    double GetCurrentStaticData(Data& rStaticData, const ConstitutiveInputMap& rConstitutiveInput) const;

    // parameters /////////////////////////////////////////////////////////////

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;


private:
    double mEnduranceStress = 0.;
    double mFatigueParameter = 0.;
};

} // namespace NuTo
