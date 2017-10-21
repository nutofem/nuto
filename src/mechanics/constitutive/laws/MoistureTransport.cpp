#include "mechanics/constitutive/laws/MoistureTransport.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"

#include <limits>

template <int TDim>
void NuTo::MoistureTransport::Evaluate(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                       const NuTo::ConstitutiveOutputMap& rConstitutiveOutput, Data& rStaticData)
{
    auto& staticData = rStaticData.GetData();
    // Copy input data to input struct
    InputData<TDim> inputData;
    for (auto& itInput : rConstitutiveInput)
    {
        switch (itInput.first)
        {
        case NuTo::Constitutive::eInput::RELATIVE_HUMIDITY:
            inputData.mRelativeHumidity = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::eInput::RELATIVE_HUMIDITY_DT1:
            inputData.mRelativeHumidity_dt1 = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::eInput::RELATIVE_HUMIDITY_GRADIENT:
            inputData.mRelativeHumidity_Gradient =
                    static_cast<ConstitutiveVector<TDim>*>(itInput.second.get())->AsVector();
            break;

        case NuTo::Constitutive::eInput::WATER_VOLUME_FRACTION:
            inputData.mWaterVolumeFraction = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::eInput::WATER_VOLUME_FRACTION_DT1:
            inputData.mWaterVolumeFraction_dt1 = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::eInput::WATER_VOLUME_FRACTION_GRADIENT:
            inputData.mWaterVolumeFraction_Gradient =
                    static_cast<ConstitutiveVector<TDim>*>(itInput.second.get())->AsVector();
            break;

        case NuTo::Constitutive::eInput::CALCULATE_STATIC_DATA:
        case NuTo::Constitutive::eInput::TIME_STEP:
            break;

        default:
            continue;
        }
    }


    // evaluate outputs
    for (const auto& itOutput : rConstitutiveOutput)
    {

        switch (itOutput.first)
        {

        // ----------------------------------------------------------------------------------------
        // Internal Gradient
        //-----------------------------------------------------------------------------------------

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
        {
            // Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first, __PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mRelativeHumidity_Gradient,
                                                    std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, TDim, 1>& internalGradientRH_B =
                    (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get())).AsVector();
            internalGradientRH_B =
                    mDiffusionCoefficientRH *
                    std::pow(1 - (inputData.mWaterVolumeFraction / mPoreVolumeFraction), mDiffusionExponentRH) *
                    inputData.mRelativeHumidity_Gradient;
        }
        break;


        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity != std::numeric_limits<double>::min());
            assert(inputData.mRelativeHumidity_dt1 != std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction_dt1 != std::numeric_limits<double>::min());

            // Calculation
            auto currentSorptionCoeff = staticData.GetCurrentSorptionCoeff();
            Eigen::Matrix<double, 1, 1>& internalGradientRH_N =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientRH_N(0, 0) =
                    mDensitySaturatedWaterVapor * (mPoreVolumeFraction - inputData.mWaterVolumeFraction) *
                            inputData.mRelativeHumidity_dt1 -
                    mDensitySaturatedWaterVapor * inputData.mWaterVolumeFraction_dt1 * inputData.mRelativeHumidity -
                    mMassExchangeRate * inputData.mWaterVolumeFraction +
                    mMassExchangeRate * inputData.mRelativeHumidity *
                            (currentSorptionCoeff(0) + currentSorptionCoeff(1) * inputData.mRelativeHumidity +
                             currentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
        }
        break;

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
        {
            // Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first, __PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mWaterVolumeFraction_Gradient,
                                                    std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, TDim, 1>& internalGradientWV_B =
                    (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get())).AsVector();
            internalGradientWV_B = mDiffusionCoefficientWV * inputData.mWaterVolumeFraction_Gradient *
                                   pow(inputData.mWaterVolumeFraction / mPoreVolumeFraction, mDiffusionExponentWV);
        }
        break;

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity != std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction_dt1 != std::numeric_limits<double>::min());

            // Calculation
            auto currentSorptionCoeff = staticData.GetCurrentSorptionCoeff();
            Eigen::Matrix<double, 1, 1>& internalGradientWV_N =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientWV_N(0, 0) =
                    mDensityWater * inputData.mWaterVolumeFraction_dt1 +
                    mMassExchangeRate * inputData.mWaterVolumeFraction -
                    mMassExchangeRate * inputData.mRelativeHumidity *
                            (currentSorptionCoeff(0) + currentSorptionCoeff(1) * inputData.mRelativeHumidity +
                             currentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
        }
        break;

        // ----------------------------------------------------------------------------------------
        // Hessian 0 (Stiffness)
        //-----------------------------------------------------------------------------------------


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_BB_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientRH_dRH_BB_H0(0, 0) =
                    mDiffusionCoefficientRH *
                    pow(1 - (inputData.mWaterVolumeFraction / mPoreVolumeFraction), mDiffusionExponentRH);
        }
        break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity != std::numeric_limits<double>::min());
            assert(inputData.mRelativeHumidity_dt1 != std::numeric_limits<double>::min());

            // Calculation
            CalculateSorptionCurveCoefficients(
                    staticData,
                    inputData
                            .mRelativeHumidity); // VHIRTHAMTODO ---> find better place to calculate, maybe in inputData
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_NN_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            auto currentSorptionCoeff = staticData.GetCurrentSorptionCoeff();
            if (mEnableModifiedTangentialStiffness)
            {
                internalGradientRH_dRH_NN_H0(0, 0) =
                        mMassExchangeRate *
                        (currentSorptionCoeff(0) + currentSorptionCoeff(1) * inputData.mRelativeHumidity +
                         currentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
            }
            else
            {
                internalGradientRH_dRH_NN_H0(0, 0) =
                        mMassExchangeRate *
                                (currentSorptionCoeff(0) + currentSorptionCoeff(1) * inputData.mRelativeHumidity * 2.0 +
                                 currentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity *
                                         3.0) -
                        mDensitySaturatedWaterVapor * inputData.mRelativeHumidity_dt1;
            }
        }
        break;


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
        {
            // Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first, __PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mRelativeHumidity_Gradient,
                                                    std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, TDim, 1>& internalGradientRH_dWV_BN_H0 =
                    (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get())).AsVector();
            if (mEnableModifiedTangentialStiffness)
            {
                internalGradientRH_dWV_BN_H0.col(0).setZero();
            }
            else
            {
                internalGradientRH_dWV_BN_H0 =
                        inputData.mRelativeHumidity_Gradient * mDiffusionCoefficientRH * mDiffusionExponentRH /
                        mPoreVolumeFraction *
                        pow(1 - (inputData.mWaterVolumeFraction / mPoreVolumeFraction), mDiffusionExponentRH - 1.0);
            }
        }
        break;


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity_dt1 != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dWV_NN_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            if (mEnableModifiedTangentialStiffness)
            {
                internalGradientRH_dWV_NN_H0(0, 0) = -mMassExchangeRate;
            }
            else
            {
                internalGradientRH_dWV_NN_H0(0, 0) =
                        mDensitySaturatedWaterVapor * inputData.mRelativeHumidity_dt1 - mMassExchangeRate;
            }
        }
        break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity != std::numeric_limits<double>::min());

            // Calculation
            CalculateSorptionCurveCoefficients(
                    staticData,
                    inputData
                            .mRelativeHumidity); // VHIRTHAMTODO ---> find better place to calculate, maybe in inputData
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dRH_NN_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            auto currentSorptionCoeff = staticData.GetCurrentSorptionCoeff();
            if (mEnableModifiedTangentialStiffness)
            {
                internalGradientWV_dRH_NN_H0(0, 0) =
                        -mMassExchangeRate *
                        (currentSorptionCoeff(0) + currentSorptionCoeff(1) * inputData.mRelativeHumidity +
                         currentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
            }
            else
            {
                internalGradientWV_dRH_NN_H0(0, 0) =
                        -mMassExchangeRate *
                        (currentSorptionCoeff(0) + currentSorptionCoeff(1) * inputData.mRelativeHumidity * 2.0 +
                         currentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity * 3.0);
            }
        }
        break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_BB_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientWV_dWV_BB_H0(0, 0) =
                    mDiffusionCoefficientWV *
                    pow(inputData.mWaterVolumeFraction / mPoreVolumeFraction, mDiffusionExponentWV);
        }
        break;


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
        {
            // Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first, __PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mWaterVolumeFraction_Gradient,
                                                    std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, TDim, 1>& internalGradientWV_dWV_BN_H0 =
                    (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get())).AsVector();
            if (mEnableModifiedTangentialStiffness)
            {
                internalGradientWV_dWV_BN_H0.col(0).setZero();
            }
            else
            {
                internalGradientWV_dWV_BN_H0 =
                        inputData.mWaterVolumeFraction_Gradient * mDiffusionCoefficientWV * mDiffusionExponentWV /
                        mPoreVolumeFraction *
                        std::pow(inputData.mWaterVolumeFraction / mPoreVolumeFraction, mDiffusionExponentWV - 1.0);
            }
        }
        break;


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_NN_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientWV_dWV_NN_H0(0, 0) = mMassExchangeRate;
        }
        break;
        // ----------------------------------------------------------------------------------------
        // Hessian 1 (Damping)
        //-----------------------------------------------------------------------------------------

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mWaterVolumeFraction != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_NN_H1 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientRH_dRH_NN_H1(0, 0) =
                    mDensitySaturatedWaterVapor * (mPoreVolumeFraction - inputData.mWaterVolumeFraction);
        }
        break;


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity != std::numeric_limits<double>::min());

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dWV_NN_H1 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientRH_dWV_NN_H1(0, 0) = -mDensitySaturatedWaterVapor * inputData.mRelativeHumidity;
        }
        break;


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_NN_H1 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();
            internalGradientWV_dWV_NN_H1(0, 0) = mDensityWater;
        }
        break;

        // ----------------------------------------------------------------------------------------
        // Boundary - Internal Gradient
        //-----------------------------------------------------------------------------------------
        // VHIRTHAMTODO --- Check everything again --- asserts missing etc.

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_Boundary_N =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();

            double relativeHumidityBoundary = (*rConstitutiveInput.at(Constitutive::eInput::RELATIVE_HUMIDITY_BOUNDARY))
                    [0]; // mControlNode->Get(Node::eDof::RELATIVEHUMIDITY).at(0,0);
            internalGradientRH_Boundary_N(0, 0) =
                    mBoundaryDiffusionCoefficientRH * (inputData.mRelativeHumidity - relativeHumidityBoundary);
        }
        break;
        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_Boundary_N =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();

            double waterVolumeFractionBoundary = GetEquilibriumWaterVolumeFraction(
                    (*rConstitutiveInput.at(Constitutive::eInput::RELATIVE_HUMIDITY_BOUNDARY))[0],
                    staticData.GetCurrentSorptionCoeff());

            internalGradientWV_Boundary_N(0, 0) =
                    mBoundaryDiffusionCoefficientWV * (inputData.mWaterVolumeFraction - waterVolumeFractionBoundary);
        }
        break;
        // ----------------------------------------------------------------------------------------
        // Boundary - Hessian 0 (Stiffness)
        //-----------------------------------------------------------------------------------------

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_Boundary_NN_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();

            internalGradientRH_dRH_Boundary_NN_H0(0, 0) = mBoundaryDiffusionCoefficientRH;
        }
        break;


        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
        {
            // Asserts
            itOutput.second->AssertIsScalar(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_Boundary_NN_H0 =
                    (*static_cast<ConstitutiveScalar*>(itOutput.second.get())).AsScalar();

            internalGradientWV_dWV_Boundary_NN_H0(0, 0) = mBoundaryDiffusionCoefficientWV;
        }
        break;

        // ----------------------------------------------------------------------------------------
        // Static data
        //-----------------------------------------------------------------------------------------


        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {

            if (staticData.GetLastRelHumValue() > inputData.mRelativeHumidity)
            {
                staticData.SetDesorption(true);
            }
            if (staticData.GetLastRelHumValue() < inputData.mRelativeHumidity)
            {
                staticData.SetDesorption(false);
            }
            staticData.SetLastRelHumValue(inputData.mRelativeHumidity);
            staticData.SetLastSorptionCoeff(staticData.GetCurrentSorptionCoeff());
            staticData.SetLastJunctionPoint(staticData.GetCurrentJunctionPoint());

            continue;
        }

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}


//! @brief ... calculates the sorption Curve coefficients when the sorption direction has changed
void NuTo::MoistureTransport::CalculateSorptionCurveCoefficients(
        Constitutive::StaticData::DataMoistureTransport& staticData, double relativeHumidity)
{
    if (mEnableSorptionHysteresis)
    {
        double lastRelHum = staticData.GetLastRelHumValue();
        auto lastSorptionCoeff = staticData.GetLastSorptionCoeff();

        Eigen::Vector4d ITDofs;
        Eigen::Vector4d ITConstants;
        Eigen::Vector4d ITRhs;
        Eigen::Vector4d ITDelta;
        Eigen::Matrix4d Jacobi;

        // Initial Coeffs = Previous Coeffs
        staticData.SetCurrentSorptionCoeff(staticData.GetLastSorptionCoeff());
        staticData.SetCurrentJunctionPoint(staticData.GetLastJunctionPoint());

        // Getting to the junction point
        // #############################

        // adsorption
        if (!staticData.IsDesorption() and relativeHumidity > staticData.GetLastJunctionPoint())
        {
            staticData.SetCurrentSorptionCoeff(mAdsorptionCoeff);
            staticData.SetCurrentJunctionPoint(1.0);
        }

        // desorption
        if (staticData.IsDesorption() and relativeHumidity < staticData.GetLastJunctionPoint())
        {
            staticData.SetCurrentSorptionCoeff(mDesorptionCoeff);
            staticData.SetCurrentJunctionPoint(0.0);
        }

        // desorption to adsorption
        if (lastRelHum < relativeHumidity and staticData.IsDesorption())
        {
            // setting initial values
            ITDofs(0) = (lastSorptionCoeff(2) + mAdsorptionCoeff(2)) / 2.0;
            ITDofs(1) = (lastSorptionCoeff(1) + mAdsorptionCoeff(1)) / 2.0;
            ITDofs(2) = (lastSorptionCoeff(0) + mAdsorptionCoeff(0)) / 2.0;
            ITDofs(3) = lastRelHum + (1.0 - lastRelHum) / 2.0;

            // calculating the constant terms
            ITConstants(0) =
                    mGradientCorrDesorptionAdsorption * (3 * lastSorptionCoeff(2) * lastRelHum * lastRelHum +
                                                         2 * lastSorptionCoeff(1) * lastRelHum + lastSorptionCoeff(0));

            ITConstants(1) = lastSorptionCoeff(2) * lastRelHum * lastRelHum + lastSorptionCoeff(1) * lastRelHum +
                             lastSorptionCoeff(0);

            ITConstants(2) = mAdsorptionCoeff(0);

            ITConstants(3) = mAdsorptionCoeff(0);

            // Newton-Raphson-Method
            const int maxIteration = 250;
            const double tolerance = 1e-8;

            int iteration = 0;
            double residual = 1;

            while (iteration < maxIteration && residual > tolerance)
            {
                Jacobi(0, 0) = 3.0 * lastRelHum * lastRelHum;
                Jacobi(0, 1) = 2.0 * lastRelHum;
                Jacobi(0, 2) = 1.0;
                Jacobi(0, 3) = 0.0;

                Jacobi(1, 0) = lastRelHum * lastRelHum;
                Jacobi(1, 1) = lastRelHum;
                Jacobi(1, 2) = 1.0;
                Jacobi(1, 3) = 0.0;

                Jacobi(2, 0) = 3.0 * ITDofs(3) * ITDofs(3);
                Jacobi(2, 1) = 2.0 * ITDofs(3);
                Jacobi(2, 2) = 1.0;
                Jacobi(2, 3) =
                        6.0 * ITDofs(3) * (ITDofs(0) - mAdsorptionCoeff(2)) + 2.0 * (ITDofs(1) - mAdsorptionCoeff(1));

                Jacobi(3, 0) = ITDofs(3) * ITDofs(3);
                Jacobi(3, 1) = ITDofs(3);
                Jacobi(3, 2) = 1.0;
                Jacobi(3, 3) = 2.0 * ITDofs(3) * (ITDofs(0) - mAdsorptionCoeff(2)) + (ITDofs(1) - mAdsorptionCoeff(1));


                ITRhs(0) = 3.0 * ITDofs(0) * lastRelHum * lastRelHum + 2.0 * ITDofs(1) * lastRelHum + ITDofs(2) -
                           ITConstants(0);

                ITRhs(1) = ITDofs(0) * lastRelHum * lastRelHum + ITDofs(1) * lastRelHum + ITDofs(2) - ITConstants(1);

                ITRhs(2) = 3.0 * (ITDofs(0) - mAdsorptionCoeff(2)) * ITDofs(3) * ITDofs(3) +
                           2.0 * (ITDofs(1) - mAdsorptionCoeff(1)) * ITDofs(3) + ITDofs(2) - ITConstants(2);

                ITRhs(3) = (ITDofs(0) - mAdsorptionCoeff(2)) * ITDofs(3) * ITDofs(3) +
                           (ITDofs(1) - mAdsorptionCoeff(1)) * ITDofs(3) + ITDofs(2) - ITConstants(3);

                ITDelta = -Jacobi.inverse() * ITRhs;
                ITDofs += ITDelta;
                residual = ITRhs.cwiseAbs().maxCoeff();
                iteration++;
            }

            Eigen::VectorXd newSorptionCoeff;
            newSorptionCoeff << ITDofs(2), ITDofs(1), ITDofs(0);
            staticData.SetCurrentSorptionCoeff(newSorptionCoeff);
            staticData.SetCurrentJunctionPoint(ITDofs(3));
        }

        // adsorption to desorption
        if (staticData.GetLastRelHumValue() > relativeHumidity and !staticData.IsDesorption())
        {
            // setting initial values
            ITDofs(0) = (lastSorptionCoeff(2) + mDesorptionCoeff(2)) / 2.0;
            ITDofs(1) = (lastSorptionCoeff(1) + mDesorptionCoeff(1)) / 2.0;
            ITDofs(2) = (lastSorptionCoeff(0) + mDesorptionCoeff(0)) / 2.0;
            ITDofs(3) = lastRelHum + (1.0 - lastRelHum) / 2.0;

            // calculating the constant terms
            ITConstants(0) =
                    mGradientCorrAdsorptionDesorption * (3 * lastSorptionCoeff(2) * lastRelHum * lastRelHum +
                                                         2 * lastSorptionCoeff(1) * lastRelHum + lastSorptionCoeff(0));

            ITConstants(1) = lastSorptionCoeff(2) * lastRelHum * lastRelHum + lastSorptionCoeff(1) * lastRelHum +
                             lastSorptionCoeff(0);

            ITConstants(2) = mDesorptionCoeff(0);

            ITConstants(3) = mDesorptionCoeff(0);

            // Newton-Raphson-Method
            const int maxIteration = 250;
            const double tolerance = 1e-8;

            int iteration = 0;
            double residual = 1;

            while (iteration < maxIteration && residual > tolerance)
            {
                Jacobi(0, 0) = 3.0 * lastRelHum * lastRelHum;
                Jacobi(0, 1) = 2.0 * lastRelHum;
                Jacobi(0, 2) = 1.0;
                Jacobi(0, 3) = 0.0;

                Jacobi(1, 0) = lastRelHum * lastRelHum;
                Jacobi(1, 1) = lastRelHum;
                Jacobi(1, 2) = 1.0;
                Jacobi(1, 3) = 0.0;

                Jacobi(2, 0) = 3.0 * ITDofs(3) * ITDofs(3);
                Jacobi(2, 1) = 2.0 * ITDofs(3);
                Jacobi(2, 2) = 1.0;
                Jacobi(2, 3) =
                        6.0 * ITDofs(3) * (ITDofs(0) - mDesorptionCoeff(2)) + 2.0 * (ITDofs(1) - mDesorptionCoeff(1));

                Jacobi(3, 0) = ITDofs(3) * ITDofs(3);
                Jacobi(3, 1) = ITDofs(3);
                Jacobi(3, 2) = 1.0;
                Jacobi(3, 3) = 2.0 * ITDofs(3) * (ITDofs(0) - mDesorptionCoeff(2)) + (ITDofs(1) - mDesorptionCoeff(1));


                ITRhs(0) = 3.0 * ITDofs(0) * lastRelHum * lastRelHum + 2.0 * ITDofs(1) * lastRelHum + ITDofs(2) -
                           ITConstants(0);

                ITRhs(1) = ITDofs(0) * lastRelHum * lastRelHum + ITDofs(1) * lastRelHum + ITDofs(2) - ITConstants(1);

                ITRhs(2) = 3.0 * (ITDofs(0) - mDesorptionCoeff(2)) * ITDofs(3) * ITDofs(3) +
                           2.0 * (ITDofs(1) - mDesorptionCoeff(1)) * ITDofs(3) + ITDofs(2) - ITConstants(2);

                ITRhs(3) = (ITDofs(0) - mDesorptionCoeff(2)) * ITDofs(3) * ITDofs(3) +
                           (ITDofs(1) - mDesorptionCoeff(1)) * ITDofs(3) + ITDofs(2) - ITConstants(3);

                ITDelta = -Jacobi.inverse() * ITRhs;
                ITDofs += ITDelta;
                residual = ITRhs.cwiseAbs().maxCoeff();
                iteration++;
            }

            Eigen::VectorXd newSorptionCoeff;
            newSorptionCoeff << ITDofs(2), ITDofs(1), ITDofs(0);
            staticData.SetCurrentSorptionCoeff(newSorptionCoeff);
            staticData.SetCurrentJunctionPoint(ITDofs(3));

            if (staticData.GetCurrentJunctionPoint() < 0)
            {
                throw NuTo::Exception(__PRETTY_FUNCTION__, "Error calculating sorption curve - junction point < 0");
            }
        }
    }
}

void NuTo::MoistureTransport::CheckValueInLimits(std::string rCallingFunction, double rValue, double rLimLower,
                                                 double rLimUpper) const
{
    if (rValue < rLimLower || rValue > rLimUpper)
        throw NuTo::Exception(rCallingFunction, "Value(" + std::to_string(rValue) + ") exceeds limits [" +
                                                        std::to_string(rLimLower) + "," + std::to_string(rLimUpper) +
                                                        "]");
}

void NuTo::MoistureTransport::CheckValuePositive(std::string rCallingFunction, double rValue,
                                                 bool rCountZeroAsPositive) const
{
    if (rValue < 0.0 || (!rCountZeroAsPositive && rValue <= 0.0))
        throw NuTo::Exception(rCallingFunction, "Value(" + std::to_string(rValue) + ") not positive");
}

void NuTo::MoistureTransport::CheckSorptionCoefficients(std::string rCallingFunction,
                                                        Eigen::VectorXd rSorptionCoefficients) const
{
    if (rSorptionCoefficients.rows() < 3 || rSorptionCoefficients.rows() > 4)
        throw NuTo::Exception(rCallingFunction, "The vector for the desorption coefficients must have 3 or 4 "
                                                "rows. --- Polynom of 3th degree --- in case of 4 "
                                                "coefficients the constant term will be deleted");
    if (rSorptionCoefficients.rows() == 4 && rSorptionCoefficients(0) != 0.0)
        throw NuTo::Exception(rCallingFunction, "The first desorption coefficients (constant term) has to be zero");
    for (int i = 0; i < rSorptionCoefficients.rows(); ++i)
    {
        if (rSorptionCoefficients(i) != 0)
            return;
    }
    throw NuTo::Exception(rCallingFunction, "All sorption coefficients are zero!");
}


bool NuTo::MoistureTransport::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol,
                                                            int rTimeDerivative) const
{
    assert(rTimeDerivative > -1);

    switch (Node::CombineDofs(rDofRow, rDofCol))
    {
    case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
        if (rTimeDerivative < 1)
            return true;
        else
            return false;
    case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
    case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
    case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
        if (rTimeDerivative < 2)
            return true;
        else
            return false;
    default:
        return false;
    }
}


//! @brief ... check parameters of the constitutive relationship
//! if one check fails, an exception is thrwon
void NuTo::MoistureTransport::CheckParameters() const
{
    CheckAdsorptionCoefficients(mAdsorptionCoeff);
    CheckDesorptionCoefficients(mDesorptionCoeff);

    CheckBoundaryDiffusionCoefficientRH(mBoundaryDiffusionCoefficientRH);
    CheckBoundaryDiffusionCoefficientWV(mBoundaryDiffusionCoefficientWV);
    CheckDensitySaturatedWaterVapor(mDensitySaturatedWaterVapor);
    CheckDensityWater(mDensityWater);
    CheckDiffusionCoefficientRH(mDiffusionCoefficientRH);
    CheckDiffusionCoefficientWV(mDiffusionCoefficientWV);
    CheckDiffusionExponentRH(mDiffusionExponentRH);
    CheckDiffusionExponentWV(mDiffusionExponentWV);
    CheckGradientCorrDesorptionAdsorption(mGradientCorrDesorptionAdsorption);
    CheckGradientCorrAdsorptionDesorption(mGradientCorrAdsorptionDesorption);
    CheckMassExchangeRate(mMassExchangeRate);
    CheckPoreVolumeFraction(mPoreVolumeFraction);
}


NuTo::ConstitutiveInputMap
NuTo::MoistureTransport::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap& rConstitutiveOutput) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {

        // Internal Gradient
        // ----------------------------------------------------------------------------------------
        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY_GRADIENT] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY] = nullptr;
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY_DT1] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION_DT1] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION_GRADIENT] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION_DT1] = nullptr;
            break;


        // Hessian 0
        // ----------------------------------------------------------------------------------------
        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY] = nullptr;
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY_DT1] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY_GRADIENT] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY_DT1] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION_GRADIENT] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
            break;


        // Hessian 1
        // ----------------------------------------------------------------------------------------

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
            break;

        // Boundary
        // ----------------------------------------------------------------------------------------

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY] = nullptr;
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY_BOUNDARY] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
            constitutiveInputMap[Constitutive::eInput::WATER_VOLUME_FRACTION] = nullptr;
            constitutiveInputMap[Constitutive::eInput::RELATIVE_HUMIDITY_BOUNDARY] = nullptr;
            break;

        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
        case NuTo::Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
            break;

        // Misc
        // ----------------------------------------------------------------------------------------


        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            continue;
            //            ProcessUnhandledOutput(__PRETTY_FUNCTION__,itOutput.first);
            //            throw Exception(std::string("[")+__PRETTY_FUNCTION__+"] output object " +
            //            Constitutive::OutputToString(itOutput.first) + " cannot be calculated by this constitutive
            //            law.");
        }
    }

    return constitutiveInputMap;
}

NuTo::Constitutive::eConstitutiveType NuTo::MoistureTransport::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::MOISTURE_TRANSPORT;
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
bool NuTo::MoistureTransport::GetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS:
        return mEnableModifiedTangentialStiffness;

    case Constitutive::eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS:
        return mEnableSorptionHysteresis;

    default:
        throw Exception(__PRETTY_FUNCTION__, std::string("Constitutive law does not have the parameter ") +
                                                     Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

// VHIRTHAMTODO check parameters needs to be called somewhere!
//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS:
        mEnableModifiedTangentialStiffness = rValue;
        return;

    case Constitutive::eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS:
        mEnableSorptionHysteresis = rValue;
        return;

    default:
        throw Exception(__PRETTY_FUNCTION__, std::string("Constitutive law does not have the parameter ") +
                                                     Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::MoistureTransport::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH:
        return mBoundaryDiffusionCoefficientRH;

    case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV:
        return mBoundaryDiffusionCoefficientWV;

    case Constitutive::eConstitutiveParameter::DENSITY_WATER:
        return mDensityWater;

    case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH:
        return mDiffusionCoefficientRH;

    case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV:
        return mDiffusionCoefficientWV;

    case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH:
        return mDiffusionExponentRH;

    case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV:
        return mDiffusionExponentWV;

    case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION:
        return mGradientCorrAdsorptionDesorption;

    case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION:
        return mGradientCorrDesorptionAdsorption;

    case Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE:
        return mMassExchangeRate;

    case Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION:
        return mPoreVolumeFraction;

    case Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR:
        return mDensitySaturatedWaterVapor;

    default:
        throw Exception(__PRETTY_FUNCTION__, std::string("Constitutive law does not have the parameter ") +
                                                     Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH:
        CheckBoundaryDiffusionCoefficientRH(rValue);
        mBoundaryDiffusionCoefficientRH = rValue;
        return;

    case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV:
        CheckBoundaryDiffusionCoefficientWV(rValue);
        mBoundaryDiffusionCoefficientWV = rValue;
        return;

    case Constitutive::eConstitutiveParameter::DENSITY_WATER:
        CheckDensityWater(rValue);
        mDensityWater = rValue;
        return;

    case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH:
        CheckDiffusionCoefficientRH(rValue);
        mDiffusionCoefficientRH = rValue;
        return;

    case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV:
        CheckDiffusionCoefficientWV(rValue);
        mDiffusionCoefficientWV = rValue;
        return;

    case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH:
        CheckDiffusionExponentRH(rValue);
        mDiffusionExponentRH = rValue;
        return;

    case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV:
        CheckDiffusionExponentWV(rValue);
        mDiffusionExponentWV = rValue;
        return;

    case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION:
        CheckGradientCorrAdsorptionDesorption(rValue);
        mGradientCorrAdsorptionDesorption = rValue;
        return;

    case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION:
        CheckGradientCorrDesorptionAdsorption(rValue);
        mGradientCorrDesorptionAdsorption = rValue;
        return;

    case Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE:
        CheckMassExchangeRate(rValue);
        mMassExchangeRate = rValue;
        return;

    case Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION:
        CheckPoreVolumeFraction(rValue);
        mPoreVolumeFraction = rValue;
        return;

    case Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR:
        CheckDensitySaturatedWaterVapor(rValue);
        mDensitySaturatedWaterVapor = rValue;
        return;

    default:
        throw Exception(__PRETTY_FUNCTION__, std::string("Constitutive law does not have the parameter ") +
                                                     Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
Eigen::VectorXd
NuTo::MoistureTransport::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION:
        return mAdsorptionCoeff;

    case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION:
        return mDesorptionCoeff;

    default:
        throw Exception(__PRETTY_FUNCTION__, std::string("Constitutive law does not have the parameter ") +
                                                     Constitutive::ConstitutiveParameterToString(rIdentifier));
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier,
                                                           Eigen::VectorXd rValue)
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION:
    {
        CheckAdsorptionCoefficients(rValue);
        switch (rValue.rows())
        {
        case 3:
        {
            mAdsorptionCoeff = rValue;
            break;
        }
        case 4:
        {
            mAdsorptionCoeff.resize(3);
            for (int i = 0; i < 3; i++)
            {
                mAdsorptionCoeff(i) = rValue(i + 1);
            }
            break;
        }
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__, "The vector for the adsorption coefficients must have "
                                                       "3 or 4 rows. --- Polynom of 3th degree --- in case of "
                                                       "4 coefficients the constant term will be deleted");
        }
        return;
    }
    case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION:
    {
        CheckDesorptionCoefficients(rValue);
        switch (rValue.rows())
        {
        case 3:
        {
            mDesorptionCoeff = rValue;
            break;
        }
        case 4:
        {
            mDesorptionCoeff.resize(3);
            for (int i = 0; i < 3; i++)
            {
                mDesorptionCoeff(i) = rValue(i + 1);
            }
            break;
        }
        default:
        {
            throw NuTo::Exception(__PRETTY_FUNCTION__, "The vector for the desorption coefficients must have "
                                                       "3 or 4 rows. --- Polynom of 3th degree --- in case of "
                                                       "4 coefficients the constant term will be deleted");
        }
        }
        return;
    }
    default:
    {
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
    }
}


//! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
//! @param rRelativeHumidity ... relative humidity
//! @param rCoeffs ... polynomial coefficients of the sorption curve
//! @return ... equilibrium water volume fraction
double NuTo::MoistureTransport::GetEquilibriumWaterVolumeFraction(double rRelativeHumidity,
                                                                  Eigen::VectorXd rCoeffs) const
{

    if (rCoeffs.rows() < 3 || rCoeffs.rows() > 4)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "The vector for the sorption coefficients must have 3 or 4 "
                                                   "rows. --- Polynom of 3th degree --- in case of 4 "
                                                   "coefficients the constant term will be deleted");
    }
    if (rCoeffs.rows() == 3)
    {
        return rCoeffs(0) * rRelativeHumidity + rCoeffs(1) * rRelativeHumidity * rRelativeHumidity +
               rCoeffs(2) * rRelativeHumidity * rRelativeHumidity * rRelativeHumidity;
    }
    else
    {
        if (rCoeffs(0) != 0.0)
        {
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                  "The first desorption coefficients (constant term) has to be zero");
        }
        return rCoeffs(1) * rRelativeHumidity + rCoeffs(2) * rRelativeHumidity * rRelativeHumidity +
               rCoeffs(3) * rRelativeHumidity * rRelativeHumidity * rRelativeHumidity;
    }
}
