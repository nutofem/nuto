#pragma once

namespace NuTo
{

class ConstraintApplicator
{
    virtual Eigen::SparseMatrix<double> Apply(Eigen::SparseMatrix<double>& K) = 0;
    virtual Eigen::VectorXd Apply(Eigen::VectorXd& f) = 0;
    virtual Eigen::VectorXd ConvertBack(Eigen::VectorXd& u) = 0;
};

} // namespace NuTo
