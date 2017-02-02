#pragma once

namespace NuTo
{
template<int TDim>
struct EvaluateDataContinuumBoundary
{

public:
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

    //! Contact element
    Eigen::MatrixXd mMortarGapMatrix;
    Eigen::MatrixXd mMortarGapMatrixPenalty;
    Eigen::VectorXd mMortarGapVector;
    Eigen::VectorXd mJacobianbyWeight;
    Eigen::VectorXd mShapeFunctionsIntegral;
    std::unordered_map<int, int> mMappingGlobal2LocalDof;

    Eigen::MatrixXd mHessianContribution;
    Eigen::MatrixXd mForceContribution;


    // Misc
    // --------------------------------------------------------------------------------------------
    double mDetJxWeightIPxSection;
    double mDetJacobian = 0;

    double mAlpha;
};
} /* namespace NuTo */
