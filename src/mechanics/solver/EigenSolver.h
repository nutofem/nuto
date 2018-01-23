#pragma once
#include <string>
#include "Eigen/Sparse"

namespace NuTo
{

Eigen::VectorXd EigenSolver(Eigen::SparseMatrix<double> A, Eigen::VectorXd b, std::string solver);

} // namespace NuTo
