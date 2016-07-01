#pragma once

#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"

namespace NuTo
{


template<int TDim>
struct EvaluateDataContinuumBoundary
{

public:
    static constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    // Mechanics
    // --------------------------------------------------------------------------------------------

    EngineeringStrain<3> mEngineeringStrainVisualize;
    EngineeringStrain<3> mEngineeringPlasticStrainVisualize;
    EngineeringStress<3> mEngineeringStressVisualize;

    ConstitutiveScalar mDamage;

    ConstitutiveScalar mLocalEqStrain;
    ConstitutiveScalar mNonlocalParameterXi;
    ConstitutiveVector<VoigtDim> mTangentLocalEqStrainStrain;
    ConstitutiveScalar mExtrapolationError;

    BoundaryType::eType mBCType;

    // Moisture Transport
    // --------------------------------------------------------------------------------------------
    //output - internal gradient
    ConstitutiveScalar          mInternalGradientRH_Boundary_N;
    ConstitutiveScalar          mInternalGradientWV_Boundary_N;
    //output - hessian 0
    ConstitutiveScalar          mInternalGradientRH_dRH_Boundary_NN_H0;
    ConstitutiveScalar          mInternalGradientWV_dWV_Boundary_NN_H0;


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
    std::map<Node::eDof, Eigen::MatrixXd> mN;

    // Misc
    // --------------------------------------------------------------------------------------------
    double mAlpha;                  // parameter for robin boundary condition
    double mDetJxWeightIPxSection;
    double mDetJacobian = 0;

};
} /* namespace NuTo */
