//
// Created by Thomas Titscher on 10/24/16.
//
#pragma once
#include "mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
namespace NuTo
{

class GradientDamageFatigueEngineeringStress : public GradientDamageEngineeringStress
{
public:
    typedef Eigen::Vector2d StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<Eigen::Vector2d>;

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLaw<GradientDamageFatigueEngineeringStress>>(
                *this, Eigen::Vector2d(0.0, 0.0));
    }

    //! @brief Evaluate the constitutive relation in 2D
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @param rStaticData Pointer to the history data.
    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput,
                  Data& rStaticData)
    {
        auto kappaAndTangent = GetCurrentStaticData<TDim>(rStaticData, rConstitutiveInput, rConstitutiveOutput);
        GradientDamageEngineeringStress::EvaluateWithKappa<TDim>(rConstitutiveInput, rConstitutiveOutput,
                                                                 kappaAndTangent.first, kappaAndTangent.second);
    }

    //! @brief Calculates the current static data based on the given CALCULATE_STATIC_DATA input.
    //! @param rStaticData History data.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    //! @return Kappa value calculated from history data and dKappa_dNonlocalEqStrain
    template <int TDim>
    std::pair<double, double> GetCurrentStaticData(Data& rStaticData, const ConstitutiveInputMap& rConstitutiveInput,
                                                   const ConstitutiveOutputMap& rConstitutiveOutput) const;

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
    double F(double s) const
    {
        s -= mEnduranceStress;
        if (s >= 0)
            return mFatigueParameter * s / mTensileStrength;
        else
            return 0;
    }

    double dF_ds(double s) const
    {
        s -= mEnduranceStress;
        if (s >= 0)
            return mFatigueParameter / mTensileStrength;
        else
            return 0;
    }

    double k(double e, double s, const StaticDataType& rData) const
    {
        if (e > rData[0])
            return e;
        double delta_e = e - rData[1];
        return rData[0] + std::abs(delta_e) * F(s);
    }

    double dk_de(double e, double s, const StaticDataType& rData) const
    {
        if (e > rData[0])
            return 1;
        double delta_e = e - rData[1];
        if (delta_e < 0)
            return -F(s);
        else
            return F(s);
    }

    double dk_ds(double e, double s, const StaticDataType& rData) const
    {
        if (e > rData[0])
            return 0;
        double delta_e = e - rData[1];
        return std::abs(delta_e) * dF_ds(s);
    }


    double mEnduranceStress = 0.;
    double mFatigueParameter = 0.;
};

} // namespace NuTo
