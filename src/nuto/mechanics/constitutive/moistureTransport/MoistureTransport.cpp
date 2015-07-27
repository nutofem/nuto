#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidityGradient2D.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFraction.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFractionGradient2D.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

#include "nuto/mechanics/constitutive/moistureTransport/MoistureTransport.h"

#include <Math.h>

NuTo::MoistureTransport::MoistureTransport()
    :   ConstitutiveBase()
{
    SetParametersValid();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase*           NuTo::MoistureTransport::AllocateStaticDataEngineeringStress_EngineeringStrain1D    (const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataMoistureTransport();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase*           NuTo::MoistureTransport::AllocateStaticDataEngineeringStress_EngineeringStrain2D    (const NuTo::ElementBase *rElement) const
{
    return new ConstitutiveStaticDataMoistureTransport();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase*           NuTo::MoistureTransport::AllocateStaticDataEngineeringStress_EngineeringStrain3D    (const NuTo::ElementBase *rElement) const
{
    return new ConstitutiveStaticDataMoistureTransport();
}

//! @brief ... calculates the sorption Curve coefficients when the sorption direction has changed
void                                        NuTo::MoistureTransport::CalculateSorptionCurveCoefficients                         (ConstitutiveStaticDataMoistureTransport* rStaticData,
                                                                                                                                 const RelativeHumidity&     rRelativeHumidity)
{
    if(mEnableSorptionHysteresis)
    {
        double                                  LastRelHum          = rStaticData->mLastRelHumValue;
        NuTo::FullVector<double,Eigen::Dynamic> LastSorptionCoeff   = rStaticData->mLastSorptionCoeff;

        NuTo::FullVector<double,4>              ITDofs;
        NuTo::FullVector<double,4>              ITConstants;
        NuTo::FullVector<double,4>              ITRhs;
        NuTo::FullVector<double,4>              ITDelta;
        NuTo::FullMatrix<double,4,4>            Jacobi;

        // Initial Coeffs = Previous Coeffs
        rStaticData->mActualSorptionCoeff(0)     = rStaticData-> mLastSorptionCoeff(0);
        rStaticData->mActualSorptionCoeff(1)     = rStaticData-> mLastSorptionCoeff(1);
        rStaticData->mActualSorptionCoeff(2)     = rStaticData-> mLastSorptionCoeff(2);
        rStaticData->mActualJunctionPoint        = rStaticData-> mLastJunctionPoint;


        // Getting to the junction point
        // #############################

        // adsorption
        if (rStaticData->mSorptionHistoryDesorption == false && rRelativeHumidity(0) > rStaticData->mLastJunctionPoint)
        {
            rStaticData->mActualSorptionCoeff(0)     = mAdsorptionCoeff(0);
            rStaticData->mActualSorptionCoeff(1)     = mAdsorptionCoeff(1);
            rStaticData->mActualSorptionCoeff(2)     = mAdsorptionCoeff(2);

            rStaticData->mActualJunctionPoint              = 1;
        }


        // desorption
        if (rStaticData->mSorptionHistoryDesorption == true && rRelativeHumidity(0) < rStaticData->mLastJunctionPoint)
        {
            rStaticData->mActualSorptionCoeff(0)     = mDesorptionCoeff(0);
            rStaticData->mActualSorptionCoeff(1)     = mDesorptionCoeff(1);
            rStaticData->mActualSorptionCoeff(2)     = mDesorptionCoeff(2);

            rStaticData->mActualJunctionPoint              = 0;
        }


        // desorption to adsorption
        if (LastRelHum < rRelativeHumidity(0) && rStaticData->mSorptionHistoryDesorption == true)
        {
            // setting initial values
            ITDofs(0)  = (rStaticData->mLastSorptionCoeff(2) + mAdsorptionCoeff(2)) / 2.0;
            ITDofs(1)  = (rStaticData->mLastSorptionCoeff(1) + mAdsorptionCoeff(1)) / 2.0;
            ITDofs(2)  = (rStaticData->mLastSorptionCoeff(0) + mAdsorptionCoeff(0)) / 2.0;
            ITDofs(3)  = LastRelHum + (1.0 - LastRelHum) / 2.0;

            // calculating the constant terms
            ITConstants(0)  = mKa * (3 * LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                     2 * LastSorptionCoeff(1) * LastRelHum              +
                                         LastSorptionCoeff(0));

            ITConstants(1)  =            LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                         LastSorptionCoeff(1) * LastRelHum              +
                                         LastSorptionCoeff(0);

            ITConstants(2)  =            mAdsorptionCoeff(0);

            ITConstants(3)  =            mAdsorptionCoeff(0);

            // Newton-Raphson-Method

            int         maxIteration    = 250;
            double      Tolerance       = 1e-8;

            int         Iteration       = 0;
            double      Residual        = 1;

            while (Iteration < maxIteration && Residual > Tolerance)
            {
                Jacobi(0,0)             = 3 * LastRelHum * LastRelHum;
                Jacobi(0,1)             = 2 * LastRelHum;
                Jacobi(0,2)             = 1;
                Jacobi(0,3)             = 0;

                Jacobi(1,0)             =     LastRelHum * LastRelHum;
                Jacobi(1,1)             =     LastRelHum;
                Jacobi(1,2)             = 1;
                Jacobi(1,3)             = 0;

                Jacobi(2,0)             = 3 * ITDofs(3)  * ITDofs(3);
                Jacobi(2,1)             = 2 * ITDofs(3);
                Jacobi(2,2)             = 1;
                Jacobi(2,3)             = 6 * ITDofs(3) * (ITDofs(0) - mAdsorptionCoeff(2)) + 2 * (ITDofs(1) - mAdsorptionCoeff(1));

                Jacobi(3,0)             =     ITDofs(3)  * ITDofs(3);
                Jacobi(3,1)             =     ITDofs(3);
                Jacobi(3,2)             = 1;
                Jacobi(3,3)             = 2 * ITDofs(3) * (ITDofs(0) - mAdsorptionCoeff(2)) +     (ITDofs(1) - mAdsorptionCoeff(1));


                ITRhs(0)                = 3 * ITDofs(0) * LastRelHum * LastRelHum +
                                          2 * ITDofs(1) * LastRelHum +
                                              ITDofs(2) - ITConstants(0);

                ITRhs(1)                =     ITDofs(0) * LastRelHum * LastRelHum +
                                              ITDofs(1) * LastRelHum +
                                              ITDofs(2) - ITConstants(1);

                ITRhs(2)                = 3 * (ITDofs(0) - mAdsorptionCoeff(2)) * ITDofs(3)  * ITDofs(3) +
                                          2 * (ITDofs(1) - mAdsorptionCoeff(1)) * ITDofs(3) +
                                              ITDofs(2) - ITConstants(2);

                ITRhs(3)                =     (ITDofs(0) - mAdsorptionCoeff(2)) * ITDofs(3)  * ITDofs(3) +
                                              (ITDofs(1) - mAdsorptionCoeff(1)) * ITDofs(3) +
                                               ITDofs(2) - ITConstants(3);

                ITDelta = -Jacobi.Inverse() * ITRhs;
                ITDofs  += ITDelta;
                Residual = ITRhs.Abs().ColumnwiseMaxCoeff()(0,0);
                Iteration++;
            }

            rStaticData->mActualSorptionCoeff(0)    = ITDofs(2);
            rStaticData->mActualSorptionCoeff(1)    = ITDofs(1);
            rStaticData->mActualSorptionCoeff(2)    = ITDofs(0);
            rStaticData->mActualJunctionPoint       = ITDofs(3);
        }

        // adsorption to desorption
        if (rStaticData->mLastRelHumValue > rRelativeHumidity(0) && rStaticData->mSorptionHistoryDesorption == false)
        {

            // setting initial values
            ITDofs(0)  = (rStaticData->mLastSorptionCoeff(2) + mDesorptionCoeff(2)) / 2.0;
            ITDofs(1)  = (rStaticData->mLastSorptionCoeff(1) + mDesorptionCoeff(1)) / 2.0;
            ITDofs(2)  = (rStaticData->mLastSorptionCoeff(0) + mDesorptionCoeff(0)) / 2.0;
            ITDofs(3)  = LastRelHum + (1.0 - LastRelHum) / 2.0;

            // calculating the constant terms
            ITConstants(0)  = mKd * (3 * LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                     2 * LastSorptionCoeff(1) * LastRelHum              +
                                         LastSorptionCoeff(0));

            ITConstants(1)  =            LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                         LastSorptionCoeff(1) * LastRelHum              +
                                         LastSorptionCoeff(0);

            ITConstants(2)  =            mDesorptionCoeff(0);

            ITConstants(3)  =            mDesorptionCoeff(0);

            // Newton-Raphson-Method

            int         maxIteration    = 250;
            double      Tolerance       = 1e-8;

            int         Iteration       = 0;
            double      Residual        = 1;

            while (Iteration < maxIteration && Residual > Tolerance)
            {
                Jacobi(0,0)             = 3 * LastRelHum * LastRelHum;
                Jacobi(0,1)             = 2 * LastRelHum;
                Jacobi(0,2)             = 1;
                Jacobi(0,3)             = 0;

                Jacobi(1,0)             =     LastRelHum * LastRelHum;
                Jacobi(1,1)             =     LastRelHum;
                Jacobi(1,2)             = 1;
                Jacobi(1,3)             = 0;

                Jacobi(2,0)             = 3 * ITDofs(3)  * ITDofs(3);
                Jacobi(2,1)             = 2 * ITDofs(3);
                Jacobi(2,2)             = 1;
                Jacobi(2,3)             = 6 * ITDofs(3) * (ITDofs(0) - mDesorptionCoeff(2)) + 2 * (ITDofs(1) - mDesorptionCoeff(1));

                Jacobi(3,0)             =     ITDofs(3)  * ITDofs(3);
                Jacobi(3,1)             =     ITDofs(3);
                Jacobi(3,2)             = 1;
                Jacobi(3,3)             = 2 * ITDofs(3) * (ITDofs(0) - mDesorptionCoeff(2)) +     (ITDofs(1) - mDesorptionCoeff(1));


                ITRhs(0)                = 3 * ITDofs(0) * LastRelHum * LastRelHum +
                                          2 * ITDofs(1) * LastRelHum +
                                              ITDofs(2) - ITConstants(0);

                ITRhs(1)                =     ITDofs(0) * LastRelHum * LastRelHum +
                                              ITDofs(1) * LastRelHum +
                                              ITDofs(2) - ITConstants(1);

                ITRhs(2)                = 3 * (ITDofs(0) - mDesorptionCoeff(2)) * ITDofs(3)  * ITDofs(3) +
                                          2 * (ITDofs(1) - mDesorptionCoeff(1)) * ITDofs(3) +
                                               ITDofs(2) - ITConstants(2);

                ITRhs(3)                =     (ITDofs(0) - mDesorptionCoeff(2)) * ITDofs(3)  * ITDofs(3) +
                                              (ITDofs(1) - mDesorptionCoeff(1)) * ITDofs(3) +
                                               ITDofs(2) - ITConstants(3);

                ITDelta = -Jacobi.Inverse() * ITRhs;
                ITDofs  += ITDelta;
                Residual = ITRhs.Abs().ColumnwiseMaxCoeff()(0,0);
                Iteration++;
            }

            rStaticData->mActualSorptionCoeff(0)    = ITDofs(2);
            rStaticData->mActualSorptionCoeff(1)    = ITDofs(1);
            rStaticData->mActualSorptionCoeff(2)    = ITDofs(0);
            rStaticData->mActualJunctionPoint       = ITDofs(3);

            //if (rStaticData->mJunctionPoint < 0)
            //{
            //    throw NuTo::MechanicsException("[NuTo::MoistureTransport::CalculateSorptionCurveCoefficients] - Error calculating sorption curve - junction point < 0");
            //}
        }
    }
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool                                        NuTo::MoistureTransport::CheckElementCompatibility                                  (Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::ELEMENT1D:
        return true;
    case NuTo::Element::ELEMENT2D:
        return true;
    case NuTo::Element::BOUNDARYELEMENT1D:
        return true;
    case NuTo::Element::BOUNDARYMOISTURETRANSPORT1D:
        return true;
    case NuTo::Element::TRUSS1D2N:
        return true;
    default:
        return false;
    }
}

//! @brief ... check the number of adsorption coefficients
//! @param rAdsorptionCoefficients ... adsorption coefficients
void                                        NuTo::MoistureTransport::CheckAdsorptionCoefficients                                (NuTo::FullVector<double,Eigen::Dynamic> rAdsorptionCoefficients) const
{
    if (rAdsorptionCoefficients.GetNumRows() < 3 || rAdsorptionCoefficients.GetNumRows() > 4)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckAdsorptionCoefficients] The vector for the adsorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
    }
    if (rAdsorptionCoefficients.GetNumRows() == 4 && rAdsorptionCoefficients(0)!=0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckAdsorptionCoefficients] The first adsorption coefficients (constant term) has to be zero");
    }
}

//! @brief ... check the boundary surface moisture transport coefficient
//! @param ... boundary surface moisture transport coefficient
void                                        NuTo::MoistureTransport::CheckBoundarySurfaceMoistureTransportCoefficient           (double rBeta) const
{
    if (rBeta < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckBoundarySurfaceMoistureTransportCoefficient] The boundary surface moisture transport coefficient must have a non-negative value.");
    }
}


//! @brief ... check the number of desorption coefficients
//! @param rDesorptionCoefficients ... desorption coefficients
void                                        NuTo::MoistureTransport::CheckDesorptionCoefficients                                (NuTo::FullVector<double,Eigen::Dynamic> rDesorptionCoefficients) const
{
    if (rDesorptionCoefficients.GetNumRows() < 3 || rDesorptionCoefficients.GetNumRows() > 4)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckDesorptionCoefficients] The vector for the desorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
    }
    if (rDesorptionCoefficients.GetNumRows() == 4 && rDesorptionCoefficients(0)!=0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckDesorptionCoefficients] The first desorption coefficients (constant term) has to be zero");
    }
}

//! @brief ... check the gradient correction when changing from desorption to adsorption
//! @param ... gradient correction when changing from desorption to adsorption
void                                        NuTo::MoistureTransport::CheckKa                                                    (double rKa) const
{
    if (rKa < 0.0 || rKa > 1.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckKa] Ka must be a value between [0 - 1].");
    }
}

//! @brief ... check the gradient correction when changing from adsorption to desorption
//! @param ... gradient correction when changing from adsorption to desorption
void                                        NuTo::MoistureTransport::CheckKd                                                    (double rKd) const
{
    if (rKd < 0.0 || rKd > 1.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckKd] Kd must be a value between [0 - 1].");
    }
}

//! @brief ... check if the mass exchange rate is non-negative
//! @param rMassExchangeRate ... mass exchange rate
void                                        NuTo::MoistureTransport::CheckMassExchangeRate                                      (double rMassExchangeRate) const
{
    if (rMassExchangeRate < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckMassExchangeRate] The mass exchange rate must have a non-negative value.");
    }
}

//! @brief ... check parameters of the constitutive relationship
//! if one check fails, an exception is thrwon
void                                        NuTo::MoistureTransport::CheckParameters                                            () const
{
    CheckAdsorptionCoefficients(mAdsorptionCoeff);
    CheckBoundarySurfaceMoistureTransportCoefficient(mBetaRelHum);
    CheckBoundarySurfaceMoistureTransportCoefficient(mBetaWVFrac);
    CheckDesorptionCoefficients(mDesorptionCoeff);
    CheckKa(mKa);
    CheckKd(mKd);
    CheckMassExchangeRate(mR);
    CheckPorosity(mEpsP);
    CheckVaporPhaseDiffusionCoefficient(mDV);
    CheckVaporPhaseDiffusionExponent(mAlphaV);
    CheckVaporPhaseSaturationDensity(mRhoVS);
    CheckWaterPhaseDensity(mRhoW);
    CheckWaterPhaseDiffusionCoefficient(mDW);
    CheckWaterPhaseDiffusionExponent(mAlphaW);
}

//! @brief ... check if the porosity is a value between 0 and 1
//! @param rPorosity ... porosity
void                                        NuTo::MoistureTransport::CheckPorosity                                              (double rPorosity) const
{
    if (rPorosity <= 0 || rPorosity >= 1)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckPorosity] The porosity must have a value between 0 and 1.");
    }
}

//! @brief ... check if vapor phase diffusion coefficient is non-negative
//! @param rVaporPhaseDiffusionCoefficient ... vapor phase diffusion coefficient
void                                        NuTo::MoistureTransport::CheckVaporPhaseDiffusionCoefficient                        (double rVaporPhaseDiffusionCoefficient) const
{
    if (rVaporPhaseDiffusionCoefficient < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckVaporPhaseDiffusionCoefficient] The vapor phase diffusion coefficient must have a non-negative value.");
    }
}

//! @brief ... check if vapor phase diffusion exponent is non-negative
//! @param rVaporPhaseDiffusionExponent ... vapor phase diffusion exponent
void                                        NuTo::MoistureTransport::CheckVaporPhaseDiffusionExponent                           (double rVaporPhaseDiffusionExponent) const
{
    if (rVaporPhaseDiffusionExponent < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckVaporPhaseDiffusionExponent] The vapor phase diffusion exponent must have a non-negative value.");
    }
}

//! @brief ... check if vapor phase saturation density is positive and non-zero
//! @param rVaporPhaseSaturationDensity ... vapor phase saturation density
void                                        NuTo::MoistureTransport::CheckVaporPhaseSaturationDensity                           (double rVaporPhaseSaturationDensity) const
{
    if (rVaporPhaseSaturationDensity <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckVaporPhaseSaturationDensity] The vapor phase saturation density must be a non-negative, non-zero value.");
    }
}

//! @brief ... check if water phase density is positive and non-zero
//! @param rWaterPhaseDensity ... water phase density
void                                        NuTo::MoistureTransport::CheckWaterPhaseDensity                                     (double rWaterPhaseDensity) const
{
    if (rWaterPhaseDensity <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckWaterPhaseDensity] The water phase density must be a non-negative, non-zero value.");
    }
}

//! @brief ... check if water phase diffusion coefficient is non-negative
//! @param rWaterPhaseDiffusionCoefficient ... water phase diffusion coefficient
void                                        NuTo::MoistureTransport::CheckWaterPhaseDiffusionCoefficient                        (double rWaterPhaseDiffusionCoefficient) const
{
    if (rWaterPhaseDiffusionCoefficient < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckWaterPhaseDiffusionCoefficient] The water phase diffusion coefficient must have a non-negative value.");
    }
}

//! @brief ... check if water phase diffusion exponent is non-negative
//! @param rWaterPhaseDiffusionExponent ... water phase diffusion exponent
void                                        NuTo::MoistureTransport::CheckWaterPhaseDiffusionExponent                           (double rWaterPhaseDiffusionExponent) const
{
    if (rWaterPhaseDiffusionExponent < 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckWaterPhaseDiffusionExponent] The water phase diffusion exponent must have a non-negative value.");
    }
}

//! @brief ... evaluate the constitutive relation in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError                         NuTo::MoistureTransport::Evaluate1D                                                 (ElementBase* rElement, int rIp,
                                                                                                                                 const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                                                                 std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput)
{
    // get section information determining which input on the constitutive level should be used
    //const SectionBase* section(rElement->GetSection());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }



    if(rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] water volume fraction needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] relative humidity needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_D1)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] firt time derivative (velocity) of water volume fraction needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_D1)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] firt time derivative (velocity) of relative humidity needed to evaluate moisture transport.");
    }
    const RelativeHumidity&         relativeHumidity            (rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)              ->second->GetRelativeHumidity());
    const RelativeHumidity&         relativeHumidityD1          (rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_D1)           ->second->GetRelativeHumidity());
    const WaterVolumeFraction&      waterVolumeFraction         (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION)          ->second->GetWaterVolumeFraction());
    const WaterVolumeFraction&      waterVolumeFractionD1       (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_D1)       ->second->GetWaterVolumeFraction());

    ConstitutiveStaticDataMoistureTransport *StaticData = (rElement->GetStaticData(rIp))->AsMoistureTransport();
    //StaticData->RelHumDecreasingHistory = false;

    for (auto itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch(itOutput->first)
        {
        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_BB:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_RH_H0_BB(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_RH_D_RH_H0_BB(0,0)     = mDV * pow(1 - (waterVolumeFraction(0) / mEpsP), mAlphaV);

            D_Residual_RH_D_RH_H0_BB.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_RH_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_RH_D_RH_H0_NN(0,0) =  mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0));
            }
            else
            {
                D_Residual_RH_D_RH_H0_NN(0,0) =  mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * 2.0 +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * 3.0) -
                                                 mRhoVS * relativeHumidityD1(0);
            }

            D_Residual_RH_D_RH_H0_NN.SetSymmetry(false);
            D_Residual_RH_D_RH_H0_NN.SetConstant(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_BN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_WV_H0_BN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_RH_D_WV_H0_BN(0,0) =  0.0;
            }
            else
            {
                const RelativeHumidity& relativeHumidityGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT)->second->GetRelativeHumidity());

                D_Residual_RH_D_WV_H0_BN(0,0) =  relativeHumidityGradient(0) * mDV * mAlphaV / mEpsP * pow(1 - (waterVolumeFraction(0) / mEpsP), mAlphaV - 1.0);
            }

            D_Residual_RH_D_WV_H0_BN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_WV_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_RH_D_WV_H0_NN(0,0) =  -mR;
            }
            else
            {
                D_Residual_RH_D_WV_H0_NN(0,0) =  -mR + mRhoVS * relativeHumidityD1(0);
            }

            D_Residual_RH_D_WV_H0_NN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_RH_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_RH_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_WV_D_RH_H0_NN(0,0) = -mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0));
            }
            else
            {
                D_Residual_WV_D_RH_H0_NN(0,0) = -mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * 2.0 +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * 3.0);
            }

            D_Residual_WV_D_RH_H0_NN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BB:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_WV_H0_BB(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_WV_D_WV_H0_BB(0,0)     = mDW * pow(waterVolumeFraction(0) / mEpsP, mAlphaW);

            D_Residual_WV_D_WV_H0_BB.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_WV_H0_BN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_WV_D_WV_H0_BN(0,0) = 0.0;
            }
            else
            {
                const WaterVolumeFraction&      waterVolumeFractionGradient (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT) ->second->GetWaterVolumeFraction());

                D_Residual_WV_D_WV_H0_BN(0,0) = waterVolumeFractionGradient(0) * mDW * mAlphaW / mEpsP * pow(waterVolumeFraction(0) / mEpsP, mAlphaW - 1.0);
            }

            D_Residual_WV_D_WV_H0_BN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_WV_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_WV_D_WV_H0_NN(0,0)=   mR;

            D_Residual_WV_D_WV_H0_NN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H1_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_RH_H1_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_RH_D_RH_H1_NN(0,0)     = mRhoVS * (mEpsP  - waterVolumeFraction(0));

            D_Residual_RH_D_RH_H1_NN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H1_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_WV_H1_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_RH_D_WV_H1_NN(0,0)     = - mRhoVS * relativeHumidity(0);

            D_Residual_RH_D_WV_H1_NN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H1_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_WV_H1_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_WV_D_WV_H1_NN(0,0)     = mRhoW;

            D_Residual_WV_D_WV_H1_NN.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_N:
        {
            ConstitutiveTangentLocal<1,1>& residualWaterPhaseN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            residualWaterPhaseN(0,0)=   mRhoW * waterVolumeFractionD1(0) +
                                        mR    * waterVolumeFraction(0) -
                                        mR    * (StaticData->mActualSorptionCoeff(0) * relativeHumidity(0) +
                                                 StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * relativeHumidity(0) +
                                                 StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * relativeHumidity(0));

            residualWaterPhaseN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_B:
        {
            const WaterVolumeFraction&      waterVolumeFractionGradient (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT) ->second->GetWaterVolumeFraction());
            ConstitutiveTangentLocal<1,1>& residualWaterPhaseB(itOutput->second->AsConstitutiveTangentLocal_1x1());

            residualWaterPhaseB(0,0)= mDW * pow(waterVolumeFraction(0) / mEpsP, mAlphaW) * waterVolumeFractionGradient(0);

            residualWaterPhaseB.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_N:
        {
            ConstitutiveTangentLocal<1,1>& residualVaporPhaseN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            residualVaporPhaseN(0,0)=   mRhoVS * (mEpsP - waterVolumeFraction(0)) * relativeHumidityD1(0) -
                                        mRhoVS * relativeHumidity(0)             * waterVolumeFractionD1(0) -
                                        mR     * waterVolumeFraction(0) +
                                        mR     * (StaticData->mActualSorptionCoeff(0) * relativeHumidity(0) +
                                                  StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * relativeHumidity(0) +
                                                  StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * relativeHumidity(0));

            residualVaporPhaseN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_B:
        {
            const RelativeHumidity& relativeHumidityGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT)->second->GetRelativeHumidity());
            ConstitutiveTangentLocal<1,1>& residualVaporPhaseB(itOutput->second->AsConstitutiveTangentLocal_1x1());

            residualVaporPhaseB(0,0)= mDV * pow(1 - (waterVolumeFraction(0) / mEpsP), mAlphaV) * relativeHumidityGradient(0);
            residualVaporPhaseB.SetSymmetry(false);
            break;
        }
        case NuTo::Constitutive::Output::BOUNDARY_SURFACE_WATER_PHASE_RESIDUAL:
        {

            const WaterVolumeFraction&     waterVolumeFractionBoundary    (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_BOUNDARY)->second->GetWaterVolumeFraction());
            if(rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_BOUNDARY)==rConstitutiveInput.end())
            {
                throw MechanicsException("[NuTo::MoistureTransport::Evaluate] water volume fraction of boundary surface needed to evaluate moisture transport.");
            }
            ConstitutiveTangentLocal<1,1>& residualBoundarySurfaceWaterPhase(itOutput->second->AsConstitutiveTangentLocal_1x1());

            residualBoundarySurfaceWaterPhase(0,0) = mBetaWVFrac*(waterVolumeFraction(0) - waterVolumeFractionBoundary(0));
            residualBoundarySurfaceWaterPhase.SetSymmetry(false);
            break;
        }
        case NuTo::Constitutive::Output::BOUNDARY_SURFACE_VAPOR_PHASE_RESIDUAL:
        {
            const RelativeHumidity&     relativeHumidityBoundary    (rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_BOUNDARY)->second->GetRelativeHumidity());
            if(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_BOUNDARY)==rConstitutiveInput.end())
            {
                throw MechanicsException("[NuTo::MoistureTransport::Evaluate] relative humidity of boundary surface needed to evaluate moisture transport.");
            }
            ConstitutiveTangentLocal<1,1>& residualBoundarySurfaceVaporPhase(itOutput->second->AsConstitutiveTangentLocal_1x1());
            residualBoundarySurfaceVaporPhase(0,0) = mBetaRelHum*(relativeHumidity(0) - relativeHumidityBoundary(0));
            residualBoundarySurfaceVaporPhase.SetSymmetry(false);
            break;
        }
        case NuTo::Constitutive::Output::BOUNDARY_SURFACE_RELATIVE_HUMIDIY_TRANSPORT_COEFFICIENT:
        {
            ConstitutiveTangentLocal<1,1>& BoundarySurfaceMoistureTransportCoefficient(itOutput->second->AsConstitutiveTangentLocal_1x1());
            BoundarySurfaceMoistureTransportCoefficient(0,0) = mBetaRelHum;
            BoundarySurfaceMoistureTransportCoefficient.SetSymmetry(false);
            break;
        }
        case NuTo::Constitutive::Output::BOUNDARY_SURFACE_WATER_VOLUME_FRACTION_TRANSPORT_COEFFICIENT:
        {
            ConstitutiveTangentLocal<1,1>& BoundarySurfaceMoistureTransportCoefficient(itOutput->second->AsConstitutiveTangentLocal_1x1());
            BoundarySurfaceMoistureTransportCoefficient(0,0) = mBetaWVFrac;
            BoundarySurfaceMoistureTransportCoefficient.SetSymmetry(false);
            break;
        }
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {

            if (StaticData->mLastRelHumValue > relativeHumidity(0))
            {
                StaticData->mSorptionHistoryDesorption = true;
            }
            if (StaticData->mLastRelHumValue < relativeHumidity(0))
            {
                StaticData->mSorptionHistoryDesorption = false;
            }
            StaticData->mLastRelHumValue = relativeHumidity(0);
            StaticData->mLastSorptionCoeff = StaticData->mActualSorptionCoeff;
            StaticData->mLastJunctionPoint = StaticData->mActualJunctionPoint;

            break;
        }
        default:
            throw MechanicsException(std::string("[NuTo::MoistureTransport::Evaluate1D ] output object)") +
                                     NuTo::Constitutive::OutputToString(itOutput->first) +
                                     std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
    }

    return Error::SUCCESSFUL;
}



//! @brief ... evaluate the constitutive relation in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
NuTo::Error::eError NuTo::MoistureTransport::Evaluate2D(NuTo::ElementBase *rElement,
                                                        int rIp,
                                                        const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase *> &rConstitutiveInput,
                                                        std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase *> &rConstitutiveOutput)
{
    // get section information determining which input on the constitutive level should be used
    //const SectionBase* section(rElement->GetSection());

    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
    }



    if(rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] water volume fraction needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] relative humidity needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_D1)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] firt time derivative (velocity) of water volume fraction needed to evaluate moisture transport.");
    }
    if(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_D1)==rConstitutiveInput.end())
    {
        throw MechanicsException("[NuTo::MoistureTransport::Evaluate] firt time derivative (velocity) of relative humidity needed to evaluate moisture transport.");
    }
    const RelativeHumidity&         relativeHumidity            (rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY)              ->second->GetRelativeHumidity());
    const RelativeHumidity&         relativeHumidityD1          (rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_D1)           ->second->GetRelativeHumidity());
    const WaterVolumeFraction&      waterVolumeFraction         (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION)          ->second->GetWaterVolumeFraction());
    const WaterVolumeFraction&      waterVolumeFractionD1       (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_D1)       ->second->GetWaterVolumeFraction());

    ConstitutiveStaticDataMoistureTransport *StaticData = (rElement->GetStaticData(rIp))->AsMoistureTransport();
    //StaticData->RelHumDecreasingHistory = false;


    for (auto itOutput = rConstitutiveOutput.begin(); itOutput != rConstitutiveOutput.end(); itOutput++)
    {
        switch(itOutput->first)
        {


        /*--------------------------------------*\
        |               STIFFNESS                |
        \*--------------------------------------*/

        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_BB:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_RH_H0_BB(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_RH_D_RH_H0_BB(0,0)     = mDV * pow(1 - (waterVolumeFraction(0) / mEpsP), mAlphaV);

            D_Residual_RH_D_RH_H0_BB.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_RH_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_RH_D_RH_H0_NN(0,0) =  mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0));
            }
            else
            {
                D_Residual_RH_D_RH_H0_NN(0,0) =  mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * 2.0 +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * 3.0) -
                                                 mRhoVS * relativeHumidityD1(0);
            }

            D_Residual_RH_D_RH_H0_NN.SetSymmetry(false);
            D_Residual_RH_D_RH_H0_NN.SetConstant(false);
            break;
        }


        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_BN:
        {
            ConstitutiveTangentLocal<2,1>& D_Residual_RH_D_WV_H0_BN(itOutput->second->AsConstitutiveTangentLocal_2x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_RH_D_WV_H0_BN.setZero();
            }
            else
            {
                const RelativeHumidityGradient2D& relativeHumidityGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT)->second->GetRelativeHumidityGradient2D());

                D_Residual_RH_D_WV_H0_BN  =  relativeHumidityGradient * mDV * mAlphaV / mEpsP * pow(1 - (waterVolumeFraction(0) / mEpsP), mAlphaV - 1.0);
            }

            D_Residual_RH_D_WV_H0_BN.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_WV_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_RH_D_WV_H0_NN(0,0) =  -mR;
            }
            else
            {
                D_Residual_RH_D_WV_H0_NN(0,0) =  -mR + mRhoVS * relativeHumidityD1(0);
            }

            D_Residual_RH_D_WV_H0_NN.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_RH_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_RH_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            CalculateSorptionCurveCoefficients(StaticData, relativeHumidity);
            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_WV_D_RH_H0_NN(0,0) = -mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0));
            }
            else
            {
                D_Residual_WV_D_RH_H0_NN(0,0) = -mR * (StaticData->mActualSorptionCoeff(0) +
                                                       StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * 2.0 +
                                                       StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * 3.0);
            }

            D_Residual_WV_D_RH_H0_NN.SetSymmetry(false);
            break;
        }

        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BB:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_WV_H0_BB(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_WV_D_WV_H0_BB(0,0)     = mDW * pow(waterVolumeFraction(0) / mEpsP, mAlphaW);

            D_Residual_WV_D_WV_H0_BB.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BN:
        {
            ConstitutiveTangentLocal<2,1>& D_Residual_WV_D_WV_H0_BN(itOutput->second->AsConstitutiveTangentLocal_2x1());

            if(mEnableModifiedTangentialStiffness)
            {
                D_Residual_WV_D_WV_H0_BN.setZero();
            }
            else
            {
                const WaterVolumeFractionGradient2D&      waterVolumeFractionGradient (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT) ->second->GetWaterVolumeFractionGradient2D());

                D_Residual_WV_D_WV_H0_BN = waterVolumeFractionGradient * mDW * mAlphaW / mEpsP * pow(waterVolumeFraction(0) / mEpsP, mAlphaW - 1.0);
            }

            D_Residual_WV_D_WV_H0_BN.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_WV_H0_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_WV_D_WV_H0_NN(0,0)=   mR;

            D_Residual_WV_D_WV_H0_NN.SetSymmetry(false);
            break;
        }



        /*--------------------------------------*\
        |                DAMPING                 |
        \*--------------------------------------*/


        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H1_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_RH_H1_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_RH_D_RH_H1_NN(0,0)     = mRhoVS * (mEpsP  - waterVolumeFraction(0));

            D_Residual_RH_D_RH_H1_NN.SetSymmetry(false);
            D_Residual_RH_D_RH_H1_NN.SetConstant(false);
            break;
        }


        case NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H1_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_RH_D_WV_H1_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_RH_D_WV_H1_NN(0,0)     = - mRhoVS * relativeHumidity(0);

            D_Residual_RH_D_WV_H1_NN.SetSymmetry(false);
            D_Residual_RH_D_WV_H1_NN.SetConstant(false);
            break;
        }


        case NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H1_NN:
        {
            ConstitutiveTangentLocal<1,1>& D_Residual_WV_D_WV_H1_NN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            D_Residual_WV_D_WV_H1_NN(0,0)     = mRhoW;

            D_Residual_WV_D_WV_H1_NN.SetSymmetry(false);
            D_Residual_WV_D_WV_H1_NN.SetConstant(false);
            break;
        }



        /*--------------------------------------*\
        |            INTERNAL FORCE              |
        \*--------------------------------------*/


        case NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_N:
        {
            ConstitutiveTangentLocal<1,1>& residualWaterPhaseN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            residualWaterPhaseN(0,0)=   mRhoW * waterVolumeFractionD1(0) +
                                        mR    * waterVolumeFraction(0) -
                                        mR    * (StaticData->mActualSorptionCoeff(0) * relativeHumidity(0) +
                                                 StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * relativeHumidity(0) +
                                                 StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * relativeHumidity(0));

            residualWaterPhaseN.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_B:
        {
            const WaterVolumeFractionGradient2D&      waterVolumeFractionGradient (rConstitutiveInput.find(NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT) ->second->GetWaterVolumeFractionGradient2D());
            ConstitutiveTangentLocal<2,1>& residualWaterPhaseB(itOutput->second->AsConstitutiveTangentLocal_2x1());
            residualWaterPhaseB.setZero();


            residualWaterPhaseB.AddBlock(0,0, mDW * pow(waterVolumeFraction(0) / mEpsP, mAlphaW) * waterVolumeFractionGradient);

            residualWaterPhaseB.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_N:
        {
            ConstitutiveTangentLocal<1,1>& residualVaporPhaseN(itOutput->second->AsConstitutiveTangentLocal_1x1());

            residualVaporPhaseN(0,0)=   mRhoVS * (mEpsP - waterVolumeFraction(0)) * relativeHumidityD1(0) -
                                        mRhoVS * relativeHumidity(0)             * waterVolumeFractionD1(0) -
                                        mR     * waterVolumeFraction(0) +
                                        mR     * (StaticData->mActualSorptionCoeff(0) * relativeHumidity(0) +
                                                  StaticData->mActualSorptionCoeff(1) * relativeHumidity(0) * relativeHumidity(0) +
                                                  StaticData->mActualSorptionCoeff(2) * relativeHumidity(0) * relativeHumidity(0) * relativeHumidity(0));

            residualVaporPhaseN.SetSymmetry(false);
            break;
        }


        case NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_B:
        {
            const RelativeHumidityGradient2D& relativeHumidityGradient(rConstitutiveInput.find(NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT)->second->GetRelativeHumidityGradient2D());
            ConstitutiveTangentLocal<2,1>& residualVaporPhaseB(itOutput->second->AsConstitutiveTangentLocal_2x1());
            residualVaporPhaseB.setZero();
            residualVaporPhaseB.AddBlock(0,0, mDV * pow(1 - (waterVolumeFraction(0) / mEpsP), mAlphaV) * relativeHumidityGradient);
            residualVaporPhaseB.SetSymmetry(false);
            break;
        }




        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {

            if (StaticData->mLastRelHumValue > relativeHumidity(0))
            {
                StaticData->mSorptionHistoryDesorption = true;
            }
            if (StaticData->mLastRelHumValue < relativeHumidity(0))
            {
                StaticData->mSorptionHistoryDesorption = false;
            }
            StaticData->mLastRelHumValue = relativeHumidity(0);
            StaticData->mLastSorptionCoeff = StaticData->mActualSorptionCoeff;
            StaticData->mLastJunctionPoint = StaticData->mActualJunctionPoint;

            break;
        }


        default:
        {
//            throw MechanicsException(std::string("[NuTo::MoistureTransport::Evaluate2D ] output object)") +
//                                     NuTo::Constitutive::OutputToString(itOutput->first) +
//                                     std::string(" could not be calculated, check the allocated material law and the section behavior."));
        }
        }
    }
    return Error::SUCCESSFUL;
}



//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
bool NuTo::MoistureTransport::GetVariableBool(NuTo::Constitutive::eConstitutiveVariable rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveVariable::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS:
        {
            return mEnableModifiedTangentialStiffness;
        }
        case Constitutive::eConstitutiveVariable::ENABLE_SORPTION_HYSTERESIS:
        {
            return mEnableSorptionHysteresis;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::GetVariableBool] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetVariableBool(NuTo::Constitutive::eConstitutiveVariable rIdentifier, bool rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveVariable::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS:
        {
            mEnableModifiedTangentialStiffness = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::ENABLE_SORPTION_HYSTERESIS:
        {
            mEnableSorptionHysteresis = rValue;
            SetParametersValid();
            return;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::SetVariableBool] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
double NuTo::MoistureTransport::GetVariableDouble(NuTo::Constitutive::eConstitutiveVariable rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveVariable::BOUNDARY_TRANSPORT_CONSTANT_GAS_PHASE:
        {
            return mBetaRelHum;
        }
        case Constitutive::eConstitutiveVariable::BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE:
        {
            return mBetaWVFrac;
        }
        case Constitutive::eConstitutiveVariable::DENSITY_WATER_PHASE:
        {
            return mRhoW;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_CONSTANT_GAS_PHASE:
        {
            return mDV;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_CONSTANT_WATER_PHASE:
        {
            return mDW;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_EXPONENT_GAS_PHASE:
        {
            return mAlphaV;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_EXPONENT_WATER_PHASE:
        {
            return mAlphaW;
        }
        case Constitutive::eConstitutiveVariable::GRADIENT_CORRECTION_ADSORPTION_DESORPTION:
        {
            return mKd;
        }
        case Constitutive::eConstitutiveVariable::GRADIENT_CORRECTION_DESORPTION_ADSORPTION:
        {
            return mKa;
        }
        case Constitutive::eConstitutiveVariable::MASS_EXCHANGE_RATE:
        {
            return mR;
        }
        case Constitutive::eConstitutiveVariable::POROSITY:
        {
            return mEpsP;
        }
        case Constitutive::eConstitutiveVariable::SATURATION_DENSITY_GAS_PHASE:
        {
            return mRhoVS;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::GetVariableDouble] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetVariableDouble(NuTo::Constitutive::eConstitutiveVariable rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveVariable::BOUNDARY_TRANSPORT_CONSTANT_GAS_PHASE:
        {
            mBetaRelHum = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE:
        {
            mBetaWVFrac = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::DENSITY_WATER_PHASE:
        {
            CheckWaterPhaseDensity(rValue);
            mRhoW = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_CONSTANT_GAS_PHASE:
        {
            CheckVaporPhaseDiffusionCoefficient(rValue);
            mDV=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_CONSTANT_WATER_PHASE:
        {
            CheckWaterPhaseDiffusionCoefficient(rValue);
            mDW=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_EXPONENT_GAS_PHASE:
        {
            CheckVaporPhaseDiffusionExponent(rValue);
            mAlphaV=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::DIFFUSION_EXPONENT_WATER_PHASE:
        {
            CheckWaterPhaseDiffusionExponent(rValue);
            mAlphaW=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::GRADIENT_CORRECTION_ADSORPTION_DESORPTION:
        {
            CheckKd(rValue);
            mKd = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::GRADIENT_CORRECTION_DESORPTION_ADSORPTION:
        {
            CheckKa(rValue);
            mKa = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::MASS_EXCHANGE_RATE:
        {
            CheckMassExchangeRate(rValue);
            mR = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::POROSITY:
        {
            CheckPorosity(rValue);
            mEpsP=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::SATURATION_DENSITY_GAS_PHASE:
        {
            CheckVaporPhaseSaturationDensity(rValue);
            mRhoVS=rValue;
            SetParametersValid();
            return;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::SetVariableDouble] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... gets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::MoistureTransport::GetVariableFullVectorDouble(NuTo::Constitutive::eConstitutiveVariable rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_ADSORPTION:
        {
            return mAdsorptionCoeff;
        }
        case Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_DESORPTION:
        {
            return mDesorptionCoeff;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::GetVariableFullVectorDouble] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... sets a variable of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested variable
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetVariableFullVectorDouble(NuTo::Constitutive::eConstitutiveVariable rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_ADSORPTION:
        {
            CheckAdsorptionCoefficients(rValue);
            switch (rValue.GetNumRows())
            {
                case 3:
                {
                    mAdsorptionCoeff = rValue;
                    break;
                }
                case 4:
                {
                    mAdsorptionCoeff.Resize(3);
                    for(int i=0; i<3; i++)
                    {
                        mAdsorptionCoeff(i) = rValue(i+1);
                    }
                    break;
                }
                default:
                {
                    throw NuTo::MechanicsException("[NuTo::MoistureTransport::SetVariableFullVectorDouble] The vector for the adsorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
                    break;
                }
            }
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_DESORPTION:
        {
            CheckDesorptionCoefficients(rValue);
            switch (rValue.GetNumRows())
            {
                case 3:
                {
                    mDesorptionCoeff = rValue;
                    break;
                }
                case 4:
                {
                    mDesorptionCoeff.Resize(3);
                    for(int i=0; i<3; i++)
                    {
                        mDesorptionCoeff(i) = rValue(i+1);
                    }
                    break;
                }
                default:
                {
                    throw NuTo::MechanicsException("[NuTo::MoistureTransport::SetVariableFullVectorDouble] The vector for the desorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
                    break;
                }
            }
            SetParametersValid();
            return;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::SetVariableFullVectorDouble] Constitutive law does not have the requested variable");
        }
    }
}



//! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
//! @param rRelativeHumidity ... relative humidity
//! @param rCoeffs ... polynomial coefficients of the sorption curve
//! @return ... equilibrium water volume fraction
double                                      NuTo::MoistureTransport::GetEquilibriumWaterVolumeFraction                           (double rRelativeHumidity,
                                                                                                                                  NuTo::FullVector<double,Eigen::Dynamic> rCoeffs) const
{

    if (rCoeffs.GetNumRows() < 3 || rCoeffs.GetNumRows() > 4)
    {
        throw NuTo::MechanicsException("[NuTo::MoistureTransport::GetEquilibriumWaterVolumeFraction] The vector for the sorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
    }
    if (rCoeffs.GetNumRows() == 3)
    {
        return rCoeffs(0) * rRelativeHumidity +
               rCoeffs(1) * rRelativeHumidity * rRelativeHumidity +
               rCoeffs(2) * rRelativeHumidity * rRelativeHumidity * rRelativeHumidity;

    }
    else
    {
        if (rCoeffs(0)!=0.0)
        {
            throw NuTo::MechanicsException("[NuTo::MoistureTransport::CheckDesorptionCoefficients] The first desorption coefficients (constant term) has to be zero");
        }
        return rCoeffs(1) * rRelativeHumidity +
               rCoeffs(2) * rRelativeHumidity * rRelativeHumidity +
               rCoeffs(3) * rRelativeHumidity * rRelativeHumidity * rRelativeHumidity;
    }
}




//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType       NuTo::MoistureTransport::GetType                                                    () const
{
    return NuTo::Constitutive::MOISTURE_TRANSPORT;
}


//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool                                        NuTo::MoistureTransport::HaveTmpStaticData                                          () const
{
    return false;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void                                        NuTo::MoistureTransport::Info                                                       (unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Mass exchange rate                  : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::MASS_EXCHANGE_RATE)              << "\n";
    rLogger << "    Porosity                            : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::POROSITY)                        << "\n";
    rLogger << "    Gas phase diffusion constant        : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::DIFFUSION_CONSTANT_GAS_PHASE)    << "\n";
    rLogger << "    Gas phase diffusion exponent        : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::DIFFUSION_EXPONENT_GAS_PHASE)    << "\n";
    rLogger << "    Gas phase saturation density        : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::SATURATION_DENSITY_GAS_PHASE)    << "\n";
    rLogger << "    Water phase density                 : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::DENSITY_WATER_PHASE)             << "\n";
    rLogger << "    Water phase diffusion constant      : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::DIFFUSION_CONSTANT_WATER_PHASE)  << "\n";
    rLogger << "    Water phase diffusion exponent      : " << this->GetVariableDouble(Constitutive::eConstitutiveVariable::DIFFUSION_EXPONENT_WATER_PHASE)  << "\n";
}





