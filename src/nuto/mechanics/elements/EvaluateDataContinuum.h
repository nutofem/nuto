#pragma once

namespace NuTo
{
template<int TDim>
struct EvaluateDataContinuum
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
    std::map<Node::eDof, const Eigen::MatrixXd*> mN;
    std::map<Node::eDof, Eigen::MatrixXd> mNIGA;

    //! @brief needed to maximize the number of function the ContinuumElement and ContinuumElementIGS have in common
    //! since IGA doesn't save the N-Matrices and thus there is no pointer on it
    //! @return the N-Matrix
    const Eigen::MatrixXd* GetNMatrix(Node::eDof dof)
    {
        if(mN.size() != 0)
            return mN.at(dof);
        else
            return &(mNIGA.at(dof));
    }

    // Misc
    // --------------------------------------------------------------------------------------------
    double mDetJxWeightIPxSection;
    double mDetJacobian = 0;
    double mTotalMass = 0;
};
} /* namespace NuTo */
