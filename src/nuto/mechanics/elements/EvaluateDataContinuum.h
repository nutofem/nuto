#pragma once

#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"

namespace NuTo
{


template<int TDim>
struct EvaluateDataContinuum
{

public:
    static constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    // Mechanics
    // --------------------------------------------------------------------------------------------

    EngineeringStrain<TDim> mEngineeringStrain;
    EngineeringStress<TDim> mEngineeringStress;
    ConstitutiveMatrix<VoigtDim, VoigtDim> mTangentStressStrain;

    EngineeringStrain<3> mEngineeringStrainVisualize;
    EngineeringStrain<3> mEngineeringPlasticStrainVisualize;
    EngineeringStress<3> mEngineeringStressVisualize;

    // Damage models
    // --------------------------------------------------------------------------------------------

    ConstitutiveScalar mDamage;
    ConstitutiveScalar mLocalEqStrain;
    ConstitutiveScalar mNonlocalEqStrain;
    ConstitutiveScalar mNonlocalParameterXi;
    ConstitutiveVector<VoigtDim> mTangentStressNonlocalEqStrain;
    ConstitutiveVector<VoigtDim> mTangentLocalEqStrainStrain;
    ConstitutiveScalar mExtrapolationError;

    // Phase field models
    // --------------------------------------------------------------------------------------------

    ConstitutiveScalar mElasticEnergyDensity;
    ConstitutiveVector<VoigtDim> mTangentElasticEnergyStrain;
    ConstitutiveVector<VoigtDim> mEngineeringStressDamagedPart;

    // Shrinkage (stress based)
    // --------------------------------------------------------------------------------------------

    ConstitutiveVector<VoigtDim> mEngineeringStress_dRH;
    ConstitutiveVector<VoigtDim> mEngineeringStress_dWV;

    EngineeringStrain<3> mShrinkageStrainVisualize;
    // Moisture Transport
    // --------------------------------------------------------------------------------------------

    //input
    ConstitutiveScalar          mRelativeHumidity;
    ConstitutiveScalar          mRelativeHumidity_dt1;
    ConstitutiveVector<TDim>    mRelativeHumidity_Gradient;
    ConstitutiveScalar          mWaterVolumeFraction;
    ConstitutiveScalar          mWaterVolumeFraction_dt1;
    ConstitutiveVector<TDim>    mWaterVolumeFraction_Gradient;
    //output - internal gradient
    ConstitutiveVector<TDim>    mInternalGradientRH_B;
    ConstitutiveScalar          mInternalGradientRH_N;
    ConstitutiveVector<TDim>    mInternalGradientWV_B;
    ConstitutiveScalar          mInternalGradientWV_N;
    //output - hessian 0
    ConstitutiveScalar          mInternalGradientRH_dRH_BB_H0;
    ConstitutiveScalar          mInternalGradientRH_dRH_NN_H0;
    ConstitutiveVector<TDim>    mInternalGradientRH_dWV_BN_H0;
    ConstitutiveScalar          mInternalGradientRH_dWV_NN_H0;
    ConstitutiveScalar          mInternalGradientWV_dRH_NN_H0;
    ConstitutiveScalar          mInternalGradientWV_dWV_BB_H0;
    ConstitutiveVector<TDim>    mInternalGradientWV_dWV_BN_H0;
    ConstitutiveScalar          mInternalGradientWV_dWV_NN_H0;
    //output - hessian 1
    ConstitutiveScalar          mInternalGradientRH_dRH_NN_H1;
    ConstitutiveScalar          mInternalGradientRH_dWV_NN_H1;
    ConstitutiveScalar          mInternalGradientWV_dWV_NN_H1;

    // Heat conduction
    // ------------------------------------------------------------------------
    ConstitutiveMatrix<TDim, TDim> mTangentHeatFluxTemperatureGradient;
    ConstitutiveScalar mTangentHeatTemperature;
    ConstitutiveScalar mTemperature;
    ConstitutiveScalar mHeatChange;
    ConstitutiveScalar mTemperatureChange;
    ConstitutiveVector<TDim> mHeatFlux;
    ConstitutiveVector<TDim> mTemperatureGradient;

    // Thermal strains
    ConstitutiveVector<VoigtDim> mDStressDTemperature;
    EngineeringStrain<3> mThermalStrain;

    // Nodal Values
    // --------------------------------------------------------------------------------------------

    std::map<Node::eDof, Eigen::VectorXd> mNodalValues;
    std::map<Node::eDof, Eigen::VectorXd> mNodalValues_dt1;


    // Shape Functions
    // --------------------------------------------------------------------------------------------

    //! @todo TT: fixed size dimensions could be used here, e.g.
    //! Eigen::Matrix<double, TDim,     Eigen::Dynamic> if B represents a normal gradient or
    //! Eigen::Matrix<double, VoigtDim, Eigen::Dynamic> if B calculates the engineering strain.
    //! However, it is harder to put them into one container. Maybe shared_ptr<Eigen::MatrixBase>?
    std::map<Node::eDof, Eigen::MatrixXd> mB;
    std::map<Node::eDof, const Eigen::MatrixXd*> mN;



    // Misc
    // --------------------------------------------------------------------------------------------

    double mDetJxWeightIPxSection;
    double mDetJacobian = 0;
    double mTotalMass = 0;

};

} /* namespace NuTo */
