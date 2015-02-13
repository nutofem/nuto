#ifndef MOISTURETRANSPORT_H
#define MOISTURETRANSPORT_H

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/moistureTransport/ConstitutiveStaticDataMoistureTransport.h"

// Temporäre Includes
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{


class MoistureTransport : public ConstitutiveBase
{
public:
    MoistureTransport();


    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase*                     AllocateStaticDataEngineeringStress_EngineeringStrain1D (const ElementBase* rElement) const override;

    //! @brief ... check the number of adsorption coefficients
    //! @param rAdsorptionCoefficients ... adsorption coefficients
    void                                            CheckAdsorptionCoefficients         (NuTo::FullVector<double,Eigen::Dynamic> rAdsorptionCoefficients) const;

    //! @brief ... check the number of desorption coefficients
    //! @param rDesorptionCoefficients ... desorption coefficients
    void                                            CheckDesorptionCoefficients         (NuTo::FullVector<double,Eigen::Dynamic> rDesorptionCoefficients) const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool                                    CheckElementCompatibility           (Element::eElementType rElementType) const override;

    //! @brief ... check if the mass exchange rate is non-negative
    //! @param rMassExchangeRate ... mass exchange rate
    void                                            CheckMassExchangeRate               (double rMassExchangeRate) const;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void                                    CheckParameters                     () const override;

    //! @brief ... check if the porosity is a value between 0 and 1
    //! @param rPorosity ... porosity
    void                                            CheckPorosity                       (double rPorosity) const;

    //! @brief ... check if vapor phase diffusion coefficient is non-negative
    //! @param rVaporPhaseDiffusionCoefficient ... vapor phase diffusion coefficient
    void                                            CheckVaporPhaseDiffusionCoefficient (double rVaporPhaseDiffusionCoefficient) const;

    //! @brief ... check if vapor phase diffusion exponent is non-negative
    //! @param rVaporPhaseDiffusionExponent ... vapor phase diffusion exponent
    void                                            CheckVaporPhaseDiffusionExponent    (double rVaporPhaseDiffusionExponent) const;

    //! @brief ... check if vapor phase saturation density is positive and non-zero
    //! @param rVaporPhaseSaturationDensity ... vapor phase saturation density
    void                                            CheckVaporPhaseSaturationDensity    (double rVaporPhaseSaturationDensity) const;

    //! @brief ... check if water phase density is positive and non-zero
    //! @param rWaterPhaseDensity ... water phase density
    void                                            CheckWaterPhaseDensity              (double rWaterPhaseDensity) const;

    //! @brief ... check if water phase diffusion coefficient is non-negative
    //! @param rWaterPhaseDiffusionCoefficient ... water phase diffusion coefficient
    void                                            CheckWaterPhaseDiffusionCoefficient (double rWaterPhaseDiffusionCoefficient) const;

    //! @brief ... check if water phase diffusion exponent is non-negative
    //! @param rWaterPhaseDiffusionExponent ... water phase diffusion exponent
    void                                            CheckWaterPhaseDiffusionExponent    (double rWaterPhaseDiffusionExponent) const;

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError                     Evaluate1D                          (ElementBase* rElement, int rIp,
                                                                                     const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                     std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput);

    //! @brief ... get adsorption coefficients as vector
    //! @return ... adsorption coefficients as vector
    virtual NuTo::FullVector<double,Eigen::Dynamic> GetAdsorptionCoefficients           () const override;

    //! @brief ... get desorption coefficients as vector
    //! @return ... desorption coefficients as vector
    virtual NuTo::FullVector<double,Eigen::Dynamic> GetDesorptionCoefficients           () const override;

    //! @brief ... get the gradient correction when changing from desorption to adsorption
    //! @return ... gradient correction when changing from desorption to adsorption
    virtual double                                  GetKa                               () const override;

    //! @brief ... get the gradient correction when changing from adsorption to desorption
    //! @return ... gradient correction when changing from adsorption to desorption
    virtual double                                  GetKd                               () const override;

    //! @brief ... get mass exchange rate between water phase and vapor phase
    //! @return ... mass exchange rate
    virtual double                                  GetMassExchangeRate                 () const override;

    //! @brief ... get porosity
    //! @return ... porosity
    virtual double                                  GetPorosity                         () const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType         GetType                             () const override;

    //! @brief ... get vapor phase diffusion coefficient
    //! @return ... vapor phase diffusion coefficient
    virtual double                                  GetVaporPhaseDiffusionCoefficient   () const override;

    //! @brief ... get vapor phase diffusion exponent
    //! @return ... vapor phase diffusion exponent
    virtual double                                  GetVaporPhaseDiffusionExponent      () const override;

    //! @brief ... get vapor phase saturation density
    //! @return ... vapor phase saturation density
    virtual double                                  GetVaporPhaseSaturationDensity      () const override;

    //! @brief ... get water phase density
    //! @return ... water phase density
    virtual double                                  GetWaterPhaseDensity                () const override;

    //! @brief ... get water phase diffusion coefficient
    //! @return ... water phase diffusion coefficient
    virtual double                                  GetWaterPhaseDiffusionCoefficient   () const override;

    //! @brief ... get water phase diffusion exponent
    //! @return ... water phase diffusion exponent
    virtual double                                  GetWaterPhaseDiffusionExponent      () const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool                                    HaveTmpStaticData                   () const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void                                            Info                                (unsigned short rVerboseLevel, Logger& rLogger) const;

    //! @brief ... set adsorption coefficients as vector
    //! @param ... adsorption coefficients as vector
    virtual void                                    SetAdsorptionCoefficients           (NuTo::FullVector<double,Eigen::Dynamic> rAdsorptionCoefficients) override;

    //! @brief ... set desorption coefficients as vector
    //! @param ... desorption coefficients as vector
    virtual void                                    SetDesorptionCoefficients           (NuTo::FullVector<double,Eigen::Dynamic> rDesorptionCoefficients) override;

    //! @brief ... set the gradient correction when changing from desorption to adsorption
    //! @param ... gradient correction when changing from desorption to adsorption
    virtual void                                    SetKa                               (double rKa) override;

    //! @brief ... set the gradient correction when changing from adsorption to desorption
    //! @param ... gradient correction when changing from adsorption to desorption
    virtual void                                    SetKd                               (double rKd) override;

    //! @brief ... set mass exchange rate between water phase and vapor phase
    //! @param ... mass exchange rate
    virtual void                                    SetMassExchangeRate                 (double rMassExchangeRate) override;

    //! @brief ... set porosity
    //! @param ... porosity
    virtual void                                    SetPorosity                         (double rPorosity) override;

    //! @brief ... set vapor phase diffusion coefficient
    //! @param ... vapor phase diffusion coefficient
    virtual void                                    SetVaporPhaseDiffusionCoefficient   (double rVaporPhaseDiffusionCoefficient) override;

    //! @brief ... set vapor phase diffusion exponent
    //! @param ... vapor phase diffusion exponent
    virtual void                                    SetVaporPhaseDiffusionExponent      (double rVaporPhaseDiffusionExponent) override;

    //! @brief ... set vapor phase saturation density
    //! @param ... vapor phase saturation density
    virtual void                                    SetVaporPhaseSaturationDensity      (double rVaporPhaseSaturationDensity) override;

    //! @brief ... set water phase density
    //! @param ... water phase density
    virtual void                                    SetWaterPhaseDensity                (double rWaterPhaseDensity) override;

    //! @brief ... set water phase diffusion coefficient
    //! @param ... water phase diffusion coefficient
    virtual void                                    SetWaterPhaseDiffusionCoefficient   (double rWaterPhaseDiffusionCoefficient) override;

    //! @brief ... set water phase diffusion exponent
    //! @param ... water phase diffusion exponent
    virtual void                                    SetWaterPhaseDiffusionExponent      (double rWaterPhaseDiffusionExponent) override;


protected:

    //! @brief ... Coefficients of the adsorption curve
    FullVector<double,Eigen::Dynamic>    mAdsorptionCoeff        {{0.0, 0.0, 0.0}};

    //! @brief ... Vapor phase diffusion exponent \f$ \alpha_V \f$
    double                              mAlphaV                 = 1.0;

    //! @brief ... Water phase diffusion exponent \f$ \alpha_W \f$
    double                              mAlphaW                 = 1.0;

    //! @brief ... Coefficients of the desorption curve
    FullVector<double,Eigen::Dynamic>   mDesorptionCoeff        {{0.0, 0.0, 0.0}};

    //! @brief ... Vapor phase diffusion coefficient \f$ D_v \f$
    double                              mDV                     = 1.0;

    //! @brief ... Water phase diffusion coefficient \f$ D_w \f$
    double                              mDW                     = 1.0;

    //! @brief ... Porosity of the specimen \f$ \Epsilon_p \f$
    double                              mEpsP                   = 0.5;

    //! @brief ... gradient correction when switching to adsorption
    double                              mKa                     = 0.0;

    //! @brief ... gradient correction when switching to desorption
    double                              mKd                     = 0.0;

    //! @brief ... Mass exchange rate between vapor phase and water phase \f$ R \f$
    double                              mR                      = 1.0;

    //! @brief ... Water phase density \f$ \rho_w \f$
    double                              mRhoW                   = 1.0;

    //! @brief ... Vapor phase saturation density \f$ \rho_v \f$
    double                  mRhoVS                  = 1.0;

};


}



#endif // MOISTURETRANSPORT_H
