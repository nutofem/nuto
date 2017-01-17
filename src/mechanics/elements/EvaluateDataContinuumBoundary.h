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

    // Misc
    // --------------------------------------------------------------------------------------------
    double mDetJxWeightIPxSection;
    double mDetJacobian = 0;

    double m1DivAlpha;
};
} /* namespace NuTo */
