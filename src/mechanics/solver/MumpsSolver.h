#pragma once

#include "Eigen/SparseCore"

namespace NuTo
{
Eigen::VectorXd MumpsSolver(Eigen::SparseMatrix<double> A, Eigen::VectorXd b, std::string solver);
}
