#include "nuto/mechanics/constitutive/laws/MoistureTransport.h"


#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"

#include <limits>
#include <Math.h>







//! @brief ... evaluate the constitutive relation
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
//! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
template <int TDim>
NuTo::Error::eError NuTo::MoistureTransport::EvaluateMoistureTransport( NuTo::ElementBase *rElement,
                                                                        int rIp,
                                                                        const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                        const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{

    ConstitutiveStaticDataMoistureTransport *StaticData = (rElement->GetStaticData(rIp))->AsMoistureTransport();

    // Copy input data to input struct
    InputData<TDim> inputData;
    for (auto itInput : rConstitutiveInput)
    {
        switch(itInput.first)
        {
        case NuTo::Constitutive::Input::RELATIVE_HUMIDITY:
            inputData.mRelativeHumidity = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::Input::RELATIVE_HUMIDITY_DT1:
            inputData.mRelativeHumidity_dt1 = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT:
            inputData.mRelativeHumidity_Gradient = (*static_cast<ConstitutiveVector<TDim>*>(itInput.second)).AsVector();
            break;

        case NuTo::Constitutive::Input::WATER_VOLUME_FRACTION:
            inputData.mWaterVolumeFraction = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_DT1:
            inputData.mWaterVolumeFraction_dt1 = (*itInput.second)[0];
            break;

        case NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT:
            inputData.mWaterVolumeFraction_Gradient = (*static_cast<ConstitutiveVector<TDim>*>(itInput.second)).AsVector();
            break;

        case NuTo::Constitutive::Input::CALCULATE_STATIC_DATA:
        case NuTo::Constitutive::Input::TIME_STEP:
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__," input object " + Constitutive::InputToString(itInput.first) + " is not needed by the moisture transport law.");
        }
    }


    // evaluate outputs
    for (auto itOutput : rConstitutiveOutput)
    {

        switch(itOutput.first)
        {

        // ----------------------------------------------------------------------------------------
        // Internal Gradient
        //-----------------------------------------------------------------------------------------



        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
        {
            //Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first,__PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mRelativeHumidity_Gradient, std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, TDim, 1>&  internalGradientRH_B = (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second)).AsVector();
            internalGradientRH_B = mDiffusionCoefficientRH* pow(1 - (inputData.mWaterVolumeFraction / mPoreVolumeFraction), mDiffusionExponentRH) * inputData.mRelativeHumidity_Gradient;
        }
            break;





        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity          !=  std::numeric_limits<double>::min());
            assert(inputData.mRelativeHumidity_dt1      !=  std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction_dt1   !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_N = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientRH_N(0,0) =     mDensitySaturatedWaterVapor * (mPoreVolumeFraction - inputData.mWaterVolumeFraction)    * inputData.mRelativeHumidity_dt1
                                          - mDensitySaturatedWaterVapor *  inputData.mWaterVolumeFraction_dt1                       * inputData.mRelativeHumidity
                                          - mMassExchangeRate           *  inputData.mWaterVolumeFraction
                                          + mMassExchangeRate           *          inputData.mRelativeHumidity * (  StaticData->mCurrentSorptionCoeff(0)
                                                                                                                  + StaticData->mCurrentSorptionCoeff(1) * inputData.mRelativeHumidity
                                                                                                                  + StaticData->mCurrentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
        }
            break;





        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
        {
            //Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first,__PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mWaterVolumeFraction_Gradient, std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, TDim, 1>&  internalGradientWV_B = (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second)).AsVector();
            internalGradientWV_B =          mDiffusionCoefficientWV * inputData.mWaterVolumeFraction_Gradient
                                          * pow(inputData.mWaterVolumeFraction / mPoreVolumeFraction, mDiffusionExponentWV);
        }
            break;





        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity          !=  std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction_dt1   !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_N = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientWV_N(0,0) =     mDensityWater       * inputData.mWaterVolumeFraction_dt1
                                          + mMassExchangeRate   * inputData.mWaterVolumeFraction
                                          - mMassExchangeRate   * inputData.mRelativeHumidity   * ( StaticData->mCurrentSorptionCoeff(0)
                                                                                                  + StaticData->mCurrentSorptionCoeff(1) * inputData.mRelativeHumidity
                                                                                                  + StaticData->mCurrentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
        }
            break;



        // ----------------------------------------------------------------------------------------
        // Hessian 0 (Stiffness)
        //-----------------------------------------------------------------------------------------



        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_BB_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientRH_dRH_BB_H0(0,0) =     mDiffusionCoefficientRH
                                                  * pow(1 - (inputData.mWaterVolumeFraction / mPoreVolumeFraction), mDiffusionExponentRH);
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity          !=  std::numeric_limits<double>::min());
            assert(inputData.mRelativeHumidity_dt1      !=  std::numeric_limits<double>::min());

            //Calculation
            CalculateSorptionCurveCoefficients(StaticData, inputData.mRelativeHumidity);    //VHIRTHAMTODO ---> find better place to calculate, maybe in inputData
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_NN_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            if(mEnableModifiedTangentialStiffness)
            {
                internalGradientRH_dRH_NN_H0(0,0) =     mMassExchangeRate   * ( StaticData->mCurrentSorptionCoeff(0)
                                                                              + StaticData->mCurrentSorptionCoeff(1) * inputData.mRelativeHumidity
                                                                              + StaticData->mCurrentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
            }
            else
            {
                internalGradientRH_dRH_NN_H0(0,0) =     mMassExchangeRate   * ( StaticData->mCurrentSorptionCoeff(0)
                                                                              + StaticData->mCurrentSorptionCoeff(1) * inputData.mRelativeHumidity * 2.0
                                                                              + StaticData->mCurrentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity * 3.0)
                                                      - mDensitySaturatedWaterVapor * inputData.mRelativeHumidity_dt1;
            }
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
        {
            //Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first,__PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mRelativeHumidity_Gradient,std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, TDim, 1>& internalGradientRH_dWV_BN_H0 = (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second)).AsVector();
            if(mEnableModifiedTangentialStiffness)
            {
                internalGradientRH_dWV_BN_H0.col(0).setZero();
            }
            else
            {
                internalGradientRH_dWV_BN_H0  =     inputData.mRelativeHumidity_Gradient    * mDiffusionCoefficientRH   * mDiffusionExponentRH  / mPoreVolumeFraction
                                                  * pow(1 - (inputData.mWaterVolumeFraction / mPoreVolumeFraction), mDiffusionExponentRH - 1.0);
            }
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity_dt1      !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dWV_NN_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            if(mEnableModifiedTangentialStiffness)
            {
                internalGradientRH_dWV_NN_H0(0,0) =   - mMassExchangeRate;
            }
            else
            {
                internalGradientRH_dWV_NN_H0(0,0) =     mDensitySaturatedWaterVapor   * inputData.mRelativeHumidity_dt1
                                                      - mMassExchangeRate;
            }
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity          !=  std::numeric_limits<double>::min());

            //Calculation
            CalculateSorptionCurveCoefficients(StaticData, inputData.mRelativeHumidity);    //VHIRTHAMTODO ---> find better place to calculate, maybe in inputData
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dRH_NN_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            if(mEnableModifiedTangentialStiffness)
            {
                internalGradientWV_dRH_NN_H0(0,0) =   - mMassExchangeRate   * ( StaticData->mCurrentSorptionCoeff(0)
                                                                              + StaticData->mCurrentSorptionCoeff(1) * inputData.mRelativeHumidity
                                                                              + StaticData->mCurrentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity);
            }
            else
            {
                internalGradientWV_dRH_NN_H0(0,0) =   - mMassExchangeRate   * ( StaticData->mCurrentSorptionCoeff(0)
                                                                              + StaticData->mCurrentSorptionCoeff(1) * inputData.mRelativeHumidity * 2.0
                                                                              + StaticData->mCurrentSorptionCoeff(2) * inputData.mRelativeHumidity * inputData.mRelativeHumidity * 3.0);
            }

        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_BB_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientWV_dWV_BB_H0(0,0) =     mDiffusionCoefficientWV
                                                  * pow(inputData.mWaterVolumeFraction / mPoreVolumeFraction, mDiffusionExponentWV);
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
        {
            //Asserts
            itOutput.second->AssertIsVector<TDim>(itOutput.first,__PRETTY_FUNCTION__);
            InputData<TDim>::AssertVectorValueIsNot(inputData.mWaterVolumeFraction_Gradient,std::numeric_limits<double>::min());
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, TDim, 1>& internalGradientWV_dWV_BN_H0 = (*static_cast<ConstitutiveVector<TDim>*>(itOutput.second)).AsVector();
            if(mEnableModifiedTangentialStiffness)
            {
                internalGradientWV_dWV_BN_H0.col(0).setZero();
            }
            else
            {
                internalGradientWV_dWV_BN_H0  =     inputData.mWaterVolumeFraction_Gradient     * mDiffusionCoefficientWV   * mDiffusionExponentWV  / mPoreVolumeFraction
                                                  * pow(inputData.mWaterVolumeFraction / mPoreVolumeFraction, mDiffusionExponentWV - 1.0);
            }
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_NN_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientWV_dWV_NN_H0(0,0) =     mMassExchangeRate;
        }
            break;





        // ----------------------------------------------------------------------------------------
        // Hessian 1 (Damping)
        //-----------------------------------------------------------------------------------------


        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mWaterVolumeFraction       !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_NN_H1 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientRH_dRH_NN_H1(0,0) =     mDensitySaturatedWaterVapor     * (mPoreVolumeFraction  - inputData.mWaterVolumeFraction);
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);
            assert(inputData.mRelativeHumidity          !=  std::numeric_limits<double>::min());

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dWV_NN_H1 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientRH_dWV_NN_H1(0,0) =   - mDensitySaturatedWaterVapor * inputData.mRelativeHumidity;
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_NN_H1 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();
            internalGradientWV_dWV_NN_H1(0,0) =     mDensityWater;
        }
            break;


        // ----------------------------------------------------------------------------------------
        // Boundary - Internal Gradient
        //-----------------------------------------------------------------------------------------

        //VHIRTHAMTODO --- Check everything again --- asserts missing etc.

        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_Boundary_N = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();

            double relativeHumidityBoundary = rElement->GetBoundaryControlNode()->GetRelativeHumidity();
            internalGradientRH_Boundary_N(0,0) =    mBoundaryDiffusionCoefficientRH * (inputData.mRelativeHumidity - relativeHumidityBoundary);
        }
            break;





        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_Boundary_N = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();

            double waterVolumeFractionBoundary = GetEquilibriumWaterVolumeFraction(rElement->GetBoundaryControlNode()->GetRelativeHumidity(),
                                                                                   StaticData->GetCurrentSorptionCoeff());
            internalGradientWV_Boundary_N(0,0) =    mBoundaryDiffusionCoefficientWV * (inputData.mWaterVolumeFraction - waterVolumeFractionBoundary);
        }
            break;



        // ----------------------------------------------------------------------------------------
        // Boundary - Hessian 0 (Stiffness)
        //-----------------------------------------------------------------------------------------

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientRH_dRH_Boundary_NN_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();

            internalGradientRH_dRH_Boundary_NN_H0(0,0) =    mBoundaryDiffusionCoefficientRH;
        }
            break;





        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
        {
            //Asserts
            itOutput.second->AssertIsScalar(itOutput.first,__PRETTY_FUNCTION__);

            //Calculation
            Eigen::Matrix<double, 1, 1>& internalGradientWV_dWV_Boundary_NN_H0 = (*static_cast<ConstitutiveScalar*>(itOutput.second)).AsScalar();

            internalGradientWV_dWV_Boundary_NN_H0(0,0) =    mBoundaryDiffusionCoefficientWV;
        }
            break;




        // ----------------------------------------------------------------------------------------
        // Static data
        //-----------------------------------------------------------------------------------------



        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {

            if (StaticData->mLastRelHumValue > inputData.mRelativeHumidity)
            {
                StaticData->mSorptionHistoryDesorption = true;
            }
            if (StaticData->mLastRelHumValue < inputData.mRelativeHumidity)
            {
                StaticData->mSorptionHistoryDesorption = false;
            }
            StaticData->mLastRelHumValue = inputData.mRelativeHumidity;
            StaticData->mLastSorptionCoeff = StaticData->mCurrentSorptionCoeff;
            StaticData->mLastJunctionPoint = StaticData->mCurrentJunctionPoint;

            break;
        }


        default:
            throw MechanicsException(__PRETTY_FUNCTION__," output object " + Constitutive::OutputToString(itOutput.first) + " could not be calculated, check the allocated material law and the section behavior.");
        }
    }
    return NuTo::Error::SUCCESSFUL;
}










//! @brief ... calculates the sorption Curve coefficients when the sorption direction has changed
void                                        NuTo::MoistureTransport::CalculateSorptionCurveCoefficients                         (ConstitutiveStaticDataMoistureTransport* rStaticData,
                                                                                                                                 double rRelativeHumidity)
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
        rStaticData->mCurrentSorptionCoeff(0)     = rStaticData-> mLastSorptionCoeff(0);
        rStaticData->mCurrentSorptionCoeff(1)     = rStaticData-> mLastSorptionCoeff(1);
        rStaticData->mCurrentSorptionCoeff(2)     = rStaticData-> mLastSorptionCoeff(2);
        rStaticData->mCurrentJunctionPoint        = rStaticData-> mLastJunctionPoint;


        // Getting to the junction point
        // #############################

        // adsorption
        if (rStaticData->mSorptionHistoryDesorption == false && rRelativeHumidity > rStaticData->mLastJunctionPoint)
        {
            rStaticData->mCurrentSorptionCoeff(0)     = mAdsorptionCoeff(0);
            rStaticData->mCurrentSorptionCoeff(1)     = mAdsorptionCoeff(1);
            rStaticData->mCurrentSorptionCoeff(2)     = mAdsorptionCoeff(2);

            rStaticData->mCurrentJunctionPoint        = 1;
        }


        // desorption
        if (rStaticData->mSorptionHistoryDesorption == true && rRelativeHumidity < rStaticData->mLastJunctionPoint)
        {
            rStaticData->mCurrentSorptionCoeff(0)     = mDesorptionCoeff(0);
            rStaticData->mCurrentSorptionCoeff(1)     = mDesorptionCoeff(1);
            rStaticData->mCurrentSorptionCoeff(2)     = mDesorptionCoeff(2);

            rStaticData->mCurrentJunctionPoint        = 0;
        }


        // desorption to adsorption
        if (LastRelHum < rRelativeHumidity && rStaticData->mSorptionHistoryDesorption == true)
        {
            // setting initial values
            ITDofs(0)  = (rStaticData->mLastSorptionCoeff(2) + mAdsorptionCoeff(2)) / 2.0;
            ITDofs(1)  = (rStaticData->mLastSorptionCoeff(1) + mAdsorptionCoeff(1)) / 2.0;
            ITDofs(2)  = (rStaticData->mLastSorptionCoeff(0) + mAdsorptionCoeff(0)) / 2.0;
            ITDofs(3)  = LastRelHum + (1.0 - LastRelHum) / 2.0;

            // calculating the constant terms
            ITConstants(0)  = mGradientCorrDesorptionAdsorption * (3 * LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                                                   2 * LastSorptionCoeff(1) * LastRelHum              +
                                                                       LastSorptionCoeff(0));

            ITConstants(1)  =                                          LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                                                       LastSorptionCoeff(1) * LastRelHum              +
                                                                       LastSorptionCoeff(0);

            ITConstants(2)  = mAdsorptionCoeff(0);

            ITConstants(3)  = mAdsorptionCoeff(0);

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

            rStaticData->mCurrentSorptionCoeff(0)    = ITDofs(2);
            rStaticData->mCurrentSorptionCoeff(1)    = ITDofs(1);
            rStaticData->mCurrentSorptionCoeff(2)    = ITDofs(0);
            rStaticData->mCurrentJunctionPoint       = ITDofs(3);
        }

        // adsorption to desorption
        if (rStaticData->mLastRelHumValue > rRelativeHumidity && rStaticData->mSorptionHistoryDesorption == false)
        {

            // setting initial values
            ITDofs(0)  = (rStaticData->mLastSorptionCoeff(2) + mDesorptionCoeff(2)) / 2.0;
            ITDofs(1)  = (rStaticData->mLastSorptionCoeff(1) + mDesorptionCoeff(1)) / 2.0;
            ITDofs(2)  = (rStaticData->mLastSorptionCoeff(0) + mDesorptionCoeff(0)) / 2.0;
            ITDofs(3)  = LastRelHum + (1.0 - LastRelHum) / 2.0;

            // calculating the constant terms
            ITConstants(0)  = mGradientCorrAdsorptionDesorption * (3 * LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                                                   2 * LastSorptionCoeff(1) * LastRelHum              +
                                                                       LastSorptionCoeff(0));

            ITConstants(1)  =                                          LastSorptionCoeff(2) * LastRelHum * LastRelHum +
                                                                       LastSorptionCoeff(1) * LastRelHum              +
                                                                       LastSorptionCoeff(0);

            ITConstants(2)  = mDesorptionCoeff(0);

            ITConstants(3)  = mDesorptionCoeff(0);

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

            rStaticData->mCurrentSorptionCoeff(0)    = ITDofs(2);
            rStaticData->mCurrentSorptionCoeff(1)    = ITDofs(1);
            rStaticData->mCurrentSorptionCoeff(2)    = ITDofs(0);
            rStaticData->mCurrentJunctionPoint       = ITDofs(3);

            if (rStaticData->mCurrentJunctionPoint < 0)
            {
                throw NuTo::MechanicsException("[NuTo::MoistureTransport::CalculateSorptionCurveCoefficients] - Error calculating sorption curve - junction point < 0");
            }
        }
    }
}

////! @brief ... check compatibility between element type and type of constitutive relationship
////! @param rElementType ... element type
////! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool                                        NuTo::MoistureTransport::CheckElementCompatibility                                  (Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::CONTINUUMELEMENT:
    case NuTo::Element::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
        return true;
    default:
        return false;
    }
}


////! @brief ... checks if the constitutive law has a specific parameter
////! @param rIdentifier ... Enum to identify the requested parameter
////! @return ... true/false
//bool NuTo::MoistureTransport::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
//{
//    switch(rIdentifier)
//    {
//        case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH:
//        case Constitutive::eConstitutiveParameter::BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE:
//        case Constitutive::eConstitutiveParameter::DENSITY_WATER:
//        case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH:
//        case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV:
//        case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH:
//        case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV:
//        case Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS:
//        case Constitutive::eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS:
//        case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION:
//        case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION:
//        case Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE:
//        case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION:
//        case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION:
//        case Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION:
//        case Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR:
//        {
//            return true;
//        }
//        default:
//        {
//            return false;
//        }
//    }
//}



//! @brief ... check parameters of the constitutive relationship
//! if one check fails, an exception is thrwon
void                                        NuTo::MoistureTransport::CheckParameters() const
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


//! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
//! @param rConstitutiveOutput ... desired constitutive outputs
//! @param rInterpolationType ... interpolation type to determine additional inputs
//! @return constitutive inputs needed for the evaluation
NuTo::ConstitutiveInputMap NuTo::MoistureTransport::GetConstitutiveInputs(  const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                                                            const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {

        // Internal Gradient
        // ----------------------------------------------------------------------------------------
        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            break;

        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY_DT1];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION_DT1];
            break;

        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT];
            break;

        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION_DT1];
            break;


        // Hessian 0
        // ----------------------------------------------------------------------------------------
        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY_DT1];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY_DT1];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
            break;


        // Hessian 1
        // ----------------------------------------------------------------------------------------

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
            break;

        // Boundary
        // ----------------------------------------------------------------------------------------

        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
            constitutiveInputMap[Constitutive::Input::RELATIVE_HUMIDITY];
            break;

        case NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
            constitutiveInputMap[Constitutive::Input::WATER_VOLUME_FRACTION];
            break;

        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
        case NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
            break;

        // Misc
        // ----------------------------------------------------------------------------------------


        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
            break;
        default:
            throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] output object " + Constitutive::OutputToString(itOutput.first) + " cannot be calculated by this constitutive law.");
        }
    }

    return constitutiveInputMap;
}


////! @brief ... checks if a constitutive law has an specific output
////! @return ... true/false
//bool NuTo::MoistureTransport::CheckOutputTypeCompatibility(NuTo::Constitutive::Output::eOutput rOutputEnum) const
//{
//    switch (rOutputEnum)
//    {
//    case Constitutive::Output::BOUNDARY_SURFACE_RELATIVE_HUMIDIY_TRANSPORT_COEFFICIENT:
//    case Constitutive::Output::BOUNDARY_SURFACE_VAPOR_PHASE_RESIDUAL:
//    case Constitutive::Output::BOUNDARY_SURFACE_WATER_PHASE_RESIDUAL:
//    case Constitutive::Output::BOUNDARY_SURFACE_WATER_VOLUME_FRACTION_TRANSPORT_COEFFICIENT:
//    case Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_BB:
//    case Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_NN:
//    case Constitutive::Output::D_RESIDUAL_RH_D_RH_H1_NN:
//    case Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_BN:
//    case Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_NN:
//    case Constitutive::Output::D_RESIDUAL_RH_D_WV_H1_NN:
//    case Constitutive::Output::D_RESIDUAL_WV_D_RH_H0_NN:
//    case Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BB:
//    case Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BN:
//    case Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_NN:
//    case Constitutive::Output::D_RESIDUAL_WV_D_WV_H1_NN:
//    case Constitutive::Output::RESIDUAL_NORM_FACTOR_RELATIVE_HUMIDITY:
//    case Constitutive::Output::RESIDUAL_NORM_FACTOR_WATER_VOLUME_FRACTION:
//    case Constitutive::Output::RESIDUAL_VAPOR_PHASE_B:
//    case Constitutive::Output::RESIDUAL_VAPOR_PHASE_N:
//    case Constitutive::Output::RESIDUAL_WATER_PHASE_B:
//    case Constitutive::Output::RESIDUAL_WATER_PHASE_N:
//    case Constitutive::Output::UPDATE_STATIC_DATA:
//    {
//        return true;
//    }
//    default:
//    {
//        return false;
//    }
//    }
//}



//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
bool NuTo::MoistureTransport::GetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS:
        {
            return mEnableModifiedTangentialStiffness;
        }
        case Constitutive::eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS:
        {
            return mEnableSorptionHysteresis;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::GetParameterBool] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetParameterBool(NuTo::Constitutive::eConstitutiveParameter rIdentifier, bool rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS:
        {
            mEnableModifiedTangentialStiffness = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS:
        {
            mEnableSorptionHysteresis = rValue;
            SetParametersValid();
            return;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::SetParameterBool] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
double NuTo::MoistureTransport::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH:
        {
            return mBoundaryDiffusionCoefficientRH;
        }
        case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV:
        {
            return mBoundaryDiffusionCoefficientWV;
        }
        case Constitutive::eConstitutiveParameter::DENSITY_WATER:
        {
            return mDensityWater;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH:
        {
            return mDiffusionCoefficientRH;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV:
        {
            return mDiffusionCoefficientWV;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH:
        {
            return mDiffusionExponentRH;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV:
        {
            return mDiffusionExponentWV;
        }
        case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION:
        {
            return mGradientCorrAdsorptionDesorption;
        }
        case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION:
        {
            return mGradientCorrDesorptionAdsorption;
        }
        case Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE:
        {
            return mMassExchangeRate;
        }
        case Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION:
        {
            return mPoreVolumeFraction;
        }
        case Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR:
        {
            return mDensitySaturatedWaterVapor;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::GetParameterDouble] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH:
        {
            CheckBoundaryDiffusionCoefficientRH(rValue);
            mBoundaryDiffusionCoefficientRH = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV:
        {
            CheckBoundaryDiffusionCoefficientWV(rValue);
            mBoundaryDiffusionCoefficientWV = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::DENSITY_WATER:
        {
            CheckDensityWater(rValue);
            mDensityWater = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH:
        {
            CheckDiffusionCoefficientRH(rValue);
            mDiffusionCoefficientRH=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV:
        {
            CheckDiffusionCoefficientWV(rValue);
            mDiffusionCoefficientWV=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH:
        {
            CheckDiffusionExponentRH(rValue);
            mDiffusionExponentRH=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV:
        {
            CheckDiffusionExponentWV(rValue);
            mDiffusionExponentWV=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION:
        {
            CheckGradientCorrAdsorptionDesorption(rValue);
            mGradientCorrAdsorptionDesorption = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION:
        {
            CheckGradientCorrDesorptionAdsorption(rValue);
            mGradientCorrDesorptionAdsorption = rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE:
        {
            CheckMassExchangeRate(rValue);
            mMassExchangeRate= rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION:
        {
            CheckPoreVolumeFraction(rValue);
            mPoreVolumeFraction=rValue;
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR:
        {
            CheckDensitySaturatedWaterVapor(rValue);
            mDensitySaturatedWaterVapor=rValue;
            SetParametersValid();
            return;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::SetParameterDouble] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... gets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @return ... value of the requested variable
NuTo::FullVector<double, Eigen::Dynamic> NuTo::MoistureTransport::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION:
        {
            return mAdsorptionCoeff;
        }
        case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION:
        {
            return mDesorptionCoeff;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::GetParameterFullVectorDouble] Constitutive law does not have the requested variable");
        }
    }
}

//! @brief ... sets a parameter of the constitutive law which is selected by an enum
//! @param rIdentifier ... Enum to identify the requested parameter
//! @param rValue ... new value for requested variable
void NuTo::MoistureTransport::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic> rValue)
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION:
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
                    throw NuTo::MechanicsException("[NuTo::MoistureTransport::SetParameterFullVectorDouble] The vector for the adsorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
                    break;
                }
            }
            SetParametersValid();
            return;
        }
        case Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION:
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
                    throw NuTo::MechanicsException("[NuTo::MoistureTransport::SetParameterFullVectorDouble] The vector for the desorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
                    break;
                }
            }
            SetParametersValid();
            return;
        }
        default:
        {
            throw MechanicsException("[NuTo::MoistureTransport::SetParameterFullVectorDouble] Constitutive law does not have the requested variable");
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



////! @brief ... print information about the object
////! @param rVerboseLevel ... verbosity of the information
//void                                        NuTo::MoistureTransport::Info                                                       (unsigned short rVerboseLevel, Logger& rLogger) const
//{
//    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
//    rLogger << "    Mass exchange rate                  : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE)              << "\n";
//    rLogger << "    PORE_VOLUME_FRACTION                            : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION)                        << "\n";
//    rLogger << "    Gas phase diffusion constant        : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH)    << "\n";
//    rLogger << "    Gas phase diffusion exponent        : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH)    << "\n";
//    rLogger << "    Gas phase saturation density        : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR)    << "\n";
//    rLogger << "    Water phase density                 : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY_WATER)             << "\n";
//    rLogger << "    Water phase diffusion constant      : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV)  << "\n";
//    rLogger << "    Water phase diffusion exponent      : " << this->GetParameterDouble(Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV)  << "\n";
//}








