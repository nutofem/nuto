#pragma once

// TODO: Replace with std::array!


#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"

// VHIRTHAMTODO Check deactivated routines, if they are still needed
// VHIRTHAMTODO Rebuild CheckXYZ routines ---> CheckParameterDouble of base class
namespace NuTo
{

enum class eDof;

// VHIRTHAMTODO make doxygen/latex description
//! @brief ... moisture transport model
//! @author Volker Hirthammer, BAM
//! @date March 2016
class MoistureTransport : public ConstitutiveBase
{
public:
    typedef Constitutive::StaticData::DataMoistureTransport StaticDataType;
    using Data = typename Constitutive::StaticData::DataContainer<StaticDataType>;

    //! @brief constructor
    MoistureTransport()
        : ConstitutiveBase()
    {
    }

    //! @brief creates corresponding IPConstitutiveLaw
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPConstitutiveLaw<MoistureTransport>>(*this, StaticDataType());
    }

private:
    //! @brief ... Input struct which stores the input data in the evaluate routine
    template <int TDim>
    struct InputData
    {
        double mRelativeHumidity = std::numeric_limits<double>::min();
        double mRelativeHumidity_dt1 = std::numeric_limits<double>::min();
        Eigen::Matrix<double, TDim, 1> mRelativeHumidity_Gradient =
                Eigen::MatrixXd::Constant(TDim, 1, std::numeric_limits<double>::min());
        double mWaterVolumeFraction = std::numeric_limits<double>::min();
        double mWaterVolumeFraction_dt1 = std::numeric_limits<double>::min();
        Eigen::Matrix<double, TDim, 1> mWaterVolumeFraction_Gradient =
                Eigen::MatrixXd::Constant(TDim, 1, std::numeric_limits<double>::min());


        static void AssertVectorValueIsNot(const Eigen::Matrix<double, TDim, 1>& rVector, double rValue)
        {
            for (unsigned int i = 0; i < TDim; ++i)
            {
                assert(rVector[i] != rValue);
            }
        }
    };

public:
    //! @brief ... evaluate the constitutive relation
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    //! @param rStaticData ... static data
    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput,
                  Data& rStaticData);


    //! @brief ... calculates the sorption Curve coefficients when the sorption direction has changed
    void CalculateSorptionCurveCoefficients(Constitutive::StaticData::DataMoistureTransport& rStaticData,
                                            double rRelativeHumidity);

private:
    //! @brief Checks if a value is within certain limits
    //! @param rCallingFunction ... name of the calling function
    //! @param rValue ... value to be checked
    //! @param rLimLower ... lower limit
    //! @param rLimUpper ... upper limit
    void CheckValueInLimits(std::string rCallingFunction, double rValue, double rLimLower, double rLimUpper) const;

    //! @brief Checks if a value is greater than zero. Throws if not.
    //! @param rCallingFunction ... name of the calling function
    //! @param rValue ... value to be checked
    //! @param rCountZeroAsPositive ... sets if zero should be counted as positive value or not
    void CheckValuePositive(std::string rCallingFunction, double rValue, bool rCountZeroAsPositive = false) const;


    //! @brief ... Checks the sorption coefficients
    //! @param rCallingFunction ... name of the calling function
    //! @param rSorptionCoefficients ... Sorption coefficients
    void CheckSorptionCoefficients(std::string rCallingFunction, Eigen::VectorXd rSorptionCoefficients) const;

public:
    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
                                               int rTimeDerivative) const override;

    //! @brief ... Checks the adsorption coefficients
    //! @param rAdsorptionCoefficients ... Adsorption coefficients
    void CheckAdsorptionCoefficients(Eigen::VectorXd rAdsorptionCoefficients) const
    {
        CheckSorptionCoefficients(__PRETTY_FUNCTION__, rAdsorptionCoefficients);
    }

    //! @brief ... Checks the desorption coefficients
    //! @param rDesorptionCoefficients ... Desorption coefficients
    void CheckDesorptionCoefficients(Eigen::VectorXd rDesorptionCoefficients) const
    {
        CheckSorptionCoefficients(__PRETTY_FUNCTION__, rDesorptionCoefficients);
    }

    //! @brief ... check the boundary diffusion coefficient of the relative humidity
    //! @param rBoundaryDiffusionCoefficientRH ... boundary diffusion coefficient of the relative humidity
    void CheckBoundaryDiffusionCoefficientRH(double rBoundaryDiffusionCoefficientRH) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rBoundaryDiffusionCoefficientRH, true);
    }

    //! @brief ... check the boundary diffusion coefficient of the water volume fraction
    //! @param rBoundaryDiffusionCoefficientWV ... boundary diffusion coefficient of the water volume fraction
    void CheckBoundaryDiffusionCoefficientWV(double rBoundaryDiffusionCoefficientWV) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rBoundaryDiffusionCoefficientWV, true);
    }

    //    //! @brief ... checks if the constitutive law has a specific parameter
    //    //! @param rIdentifier ... Enum to identify the requested parameter
    //    //! @return ... true/false
    //    virtual bool                                    CheckHaveParameter(Constitutive::eConstitutiveParameter
    //    rIdentifier) const override;

    //! @brief ... check the gradient correction when changing from desorption to adsorption
    //! @param ... gradient correction when changing from desorption to adsorption
    void CheckGradientCorrDesorptionAdsorption(double rGradientCorrDesorptionAdsorption) const
    {
        CheckValueInLimits(__PRETTY_FUNCTION__, rGradientCorrDesorptionAdsorption, 0.0, 1.0);
    }

    //! @brief ... check the gradient correction when changing from adsorption to desorption
    //! @param ... gradient correction when changing from adsorption to desorption
    void CheckGradientCorrAdsorptionDesorption(double rGradientCorrAdsorptionDesorption) const
    {
        CheckValueInLimits(__PRETTY_FUNCTION__, rGradientCorrAdsorptionDesorption, 0.0, 1.0);
    }

    //! @brief ... check if the mass exchange rate is non-negative
    //! @param rMassExchangeRate ... mass exchange rate
    void CheckMassExchangeRate(double rMassExchangeRate) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rMassExchangeRate, true);
    }

    //    //! @brief ... checks if a constitutive law has an specific output
    //    //! @return ... true/false
    //    virtual bool                                    CheckOutputTypeCompatibility
    //    (Constitutive::eOutput rOutputEnum) const override;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void CheckParameters() const override;

    //! @brief ... check if the porosity is a value between 0 and 1
    //! @param rPorosity ... porosity
    void CheckPoreVolumeFraction(double rPoreVolumeFraction) const
    {
        CheckValueInLimits(__PRETTY_FUNCTION__, rPoreVolumeFraction, 0.0, 1.0);
    }

    //! @brief ... check if the relative humidity diffusion coefficient is non-negative
    //! @param rDiffusionCoefficientRH ... relative humidity diffusion coefficient
    void CheckDiffusionCoefficientRH(double rDiffusionCoefficientRH) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rDiffusionCoefficientRH);
    }

    //! @brief ... check if the water volume fraction diffusion coefficient is non-negative
    //! @param rDiffusionCoefficientWV ... water volume fraction diffusion coefficient
    void CheckDiffusionCoefficientWV(double rDiffusionCoefficientWV) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rDiffusionCoefficientWV);
    }

    //! @brief ... check if the relative humidity diffusion exponent is non-negative
    //! @param rDiffusionExponentRH ... relative humidity diffusion exponent
    void CheckDiffusionExponentRH(double rDiffusionExponentRH) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rDiffusionExponentRH);
    }

    //! @brief ... check if the water volume fraction diffusion exponent is non-negative
    //! @param rDiffusionExponentWV ... water volume fraction diffusion exponent
    void CheckDiffusionExponentWV(double rDiffusionExponentWV) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rDiffusionExponentWV);
    }


    //! @brief ... check if the C is positive and non-zero
    //! @param rDensitySaturatedWaterVapor ... the density of saturated water vapor
    void CheckDensitySaturatedWaterVapor(double rDensitySaturatedWaterVapor) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rDensitySaturatedWaterVapor);
    }

    //! @brief ... check if water phase density is positive and non-zero
    //! @param rWaterPhaseDensity ... water phase density
    void CheckDensityWater(double rDensityWater) const
    {
        CheckValuePositive(__PRETTY_FUNCTION__, rDensityWater);
    }

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual bool GetParameterBool(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterBool(Constitutive::eConstitutiveParameter rIdentifier, bool rValue) override;

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual Eigen::VectorXd
    GetParameterFullVectorDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterFullVectorDouble(Constitutive::eConstitutiveParameter rIdentifier,
                                              Eigen::VectorXd rValue) override;


    // VHIRTHAMTODO Check if static function better?
    //! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
    //! @param rRelativeHumidity ... relative humidity
    //! @param rCoeffs ... polynomial coefficients of the sorption curve
    //! @return ... equilibrium water volume fraction
    virtual double GetEquilibriumWaterVolumeFraction(double rRelativeHumidity, Eigen::VectorXd rCoeffs) const override;


    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput) const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or
    //! stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const override
    {
        return false;
    }


protected:
    //! @brief Coefficients of the adsorption curve.
    Eigen::VectorXd mAdsorptionCoeff = Eigen::VectorXd::Zero(3);

    //! @brief Coefficients of the desorption curve.
    Eigen::VectorXd mDesorptionCoeff = Eigen::VectorXd::Zero(3);

    //! @brief Controls if a modified tangential stiffness should be used during Newton iteration (less terms to
    //! calculate in Hessian_0).
    bool mEnableModifiedTangentialStiffness = false;

    //! @brief Controls if the sorption hysteresis model should be used.
    bool mEnableSorptionHysteresis = false;

    //! @brief Boundary surface relative humidity diffusion coefficient.
    double mBoundaryDiffusionCoefficientRH = 1.0;

    //! @brief Boundary surface water volume fraction diffusion coefficient.
    double mBoundaryDiffusionCoefficientWV = 1.0;

    //! @brief Relative humidity diffusion coefficient \f$ D_v \f$.
    double mDiffusionCoefficientRH = 1.0;

    //! @brief Water phase diffusion coefficient \f$ D_w \f$.
    double mDiffusionCoefficientWV = 1.0;

    //! @brief Relative humidity diffusion exponent \f$ \alpha_V \f$.
    double mDiffusionExponentRH = 1.0;

    //! @brief Water volume fraction diffusion exponent \f$ \alpha_W \f$.
    double mDiffusionExponentWV = 1.0;

    //! @brief Pore Volume Fraction of the specimen \f$ E_p \f$.
    double mPoreVolumeFraction = 0.5;

    //! @brief gradient correction when switching to adsorption.
    double mGradientCorrDesorptionAdsorption = 0.0;

    //! @brief gradient correction when switching to desorption.
    double mGradientCorrAdsorptionDesorption = 0.0;

    //! @brief Mass exchange rate between vapor phase and water phase \f$ R \f$.
    double mMassExchangeRate = 1.0;

    //! @brief Density of water \f$ \rho_w \f$.
    double mDensityWater = 1.0;

    //! @brief Density of saturated water vapor \f$ \rho_v \f$.
    double mDensitySaturatedWaterVapor = 1.0;
};
}
