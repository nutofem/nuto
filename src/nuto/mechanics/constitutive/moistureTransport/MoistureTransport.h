#ifndef MOISTURETRANSPORT_H
#define MOISTURETRANSPORT_H

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"


namespace NuTo
{

class MoistureTransport : public ConstitutiveBase
{
public:
    MoistureTransport();

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool                                CheckElementCompatibility           (Element::eElementType rElementType) const override;

    //! @brief ... check if the mass exchange rate is non-negative
    //! @param rMassExchangeRate ... mass exchange rate
    void                                        CheckMassExchangeRate               (double rMassExchangeRate) const;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void                                CheckParameters                     () const override;

    //! @brief ... check if the porosity is a value between 0 and 1
    //! @param rPorosity ... porosity
    void                                        CheckPorosity                       (double rPorosity) const;

    //! @brief ... check if vapor phase diffusion coefficient is non-negative
    //! @param rVaporPhaseDiffusionCoefficient ... vapor phase diffusion coefficient
    void                                        CheckVaporPhaseDiffusionCoefficient (double rVaporPhaseDiffusionCoefficient) const;

    //! @brief ... check if vapor phase diffusion exponent is non-negative
    //! @param rVaporPhaseDiffusionExponent ... vapor phase diffusion exponent
    void                                        CheckVaporPhaseDiffusionExponent    (double rVaporPhaseDiffusionExponent) const;

    //! @brief ... check if vapor phase saturation density is positive and non-zero
    //! @param rVaporPhaseSaturationDensity ... vapor phase saturation density
    void                                        CheckVaporPhaseSaturationDensity    (double rVaporPhaseSaturationDensity) const;

    //! @brief ... check if water phase density is positive and non-zero
    //! @param rWaterPhaseDensity ... water phase density
    void                                        CheckWaterPhaseDensity              (double rWaterPhaseDensity) const;

    //! @brief ... check if water phase diffusion coefficient is non-negative
    //! @param rWaterPhaseDiffusionCoefficient ... water phase diffusion coefficient
    void                                        CheckWaterPhaseDiffusionCoefficient (double rWaterPhaseDiffusionCoefficient) const;

    //! @brief ... check if water phase diffusion exponent is non-negative
    //! @param rWaterPhaseDiffusionExponent ... water phase diffusion exponent
    void                                        CheckWaterPhaseDiffusionExponent    (double rWaterPhaseDiffusionExponent) const;

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError                 Evaluate1D                          (ElementBase* rElement, int rIp,
                                                                                     const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                     std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput);

    //! @brief ... get mass exchange rate between water phase and vapor phase
    //! @return ... mass exchange rate
    virtual double                              GetMassExchangeRate                 () const override;

    //! @brief ... get porosity
    //! @return ... porosity
    virtual double                              GetPorosity                         () const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType     GetType                             () const override;

    //! @brief ... get vapor phase diffusion coefficient
    //! @return ... vapor phase diffusion coefficient
    virtual double                              GetVaporPhaseDiffusionCoefficient   () const override;

    //! @brief ... get vapor phase diffusion exponent
    //! @return ... vapor phase diffusion exponent
    virtual double                              GetVaporPhaseDiffusionExponent      () const override;

    //! @brief ... get vapor phase saturation density
    //! @return ... vapor phase saturation density
    virtual double                              GetVaporPhaseSaturationDensity      () const override;

    //! @brief ... get water phase density
    //! @return ... water phase density
    virtual double                              GetWaterPhaseDensity                () const override;

    //! @brief ... get water phase diffusion coefficient
    //! @return ... water phase diffusion coefficient
    virtual double                              GetWaterPhaseDiffusionCoefficient   () const override;

    //! @brief ... get water phase diffusion exponent
    //! @return ... water phase diffusion exponent
    virtual double                              GetWaterPhaseDiffusionExponent      () const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool                                HaveTmpStaticData                   () const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void                                        Info                                (unsigned short rVerboseLevel, Logger& rLogger) const;

    //! @brief ... set mass exchange rate between water phase and vapor phase
    //! @param ... mass exchange rate
    virtual void                                SetMassExchangeRate                 (double rMassExchangeRate) override;

    //! @brief ... set porosity
    //! @param ... porosity
    virtual void                                SetPorosity                         (double rPorosity) override;

    //! @brief ... set vapor phase diffusion coefficient
    //! @param ... vapor phase diffusion coefficient
    virtual void                                SetVaporPhaseDiffusionCoefficient   (double rVaporPhaseDiffusionCoefficient) override;

    //! @brief ... set vapor phase diffusion exponent
    //! @param ... vapor phase diffusion exponent
    virtual void                                SetVaporPhaseDiffusionExponent      (double rVaporPhaseDiffusionExponent) override;

    //! @brief ... set vapor phase saturation density
    //! @param ... vapor phase saturation density
    virtual void                                SetVaporPhaseSaturationDensity      (double rVaporPhaseSaturationDensity) override;

    //! @brief ... set water phase density
    //! @param ... water phase density
    virtual void                                SetWaterPhaseDensity                (double rWaterPhaseDensity) override;

    //! @brief ... set water phase diffusion coefficient
    //! @param ... water phase diffusion coefficient
    virtual void                                SetWaterPhaseDiffusionCoefficient   (double rWaterPhaseDiffusionCoefficient) override;

    //! @brief ... set water phase diffusion exponent
    //! @param ... water phase diffusion exponent
    virtual void                                SetWaterPhaseDiffusionExponent      (double rWaterPhaseDiffusionExponent) override;


protected:

    FullVector<double,3> mActualSorptionCoeff   {{0.20, 0.0, 0.0}};

    //! @brief ... Vapor phase diffusion exponent \f$ \alpha_V \f$
    double mAlphaV  = 1.0;

    //! @brief ... Water phase diffusion exponent \f$ \alpha_W \f$
    double mAlphaW  = 1.0;

    //! @brief ... Vapor phase diffusion coefficient \f$ D_v \f$
    double mDV      = 1.0;

    //! @brief ... Water phase diffusion coefficient \f$ D_w \f$
    double mDW      = 1.0;

    //! @brief ... Porosity of the specimen \f$ \Epsilon_p \f$
    double mEpsP    = 0.5;

    //! @brief ... Mass exchange rate between vapor phase and water phase \f$ R \f$
    double mR       = 1.0;

    //! @brief ... Water phase density \f$ \rho_w \f$
    double mRhoW    = 1.0;

    //! @brief ... Vapor phase saturation density \f$ \rho_v \f$
    double mRhoVS   = 1.0;

};


}



#endif // MOISTURETRANSPORT_H
