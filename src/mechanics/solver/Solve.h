#pragma once
#include <Eigen/Core>
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"
#include "mechanics/constraints/Constraints.h"

namespace NuTo
{

Eigen::VectorXd Solve(Eigen::SparseMatrix<double> A, Eigen::VectorXd b);


// TODO: DofType, numIndependentDofs and time need to go eventually
GlobalDofVector Solve(GlobalDofMatrixSparse K, GlobalDofVector f, Constraint::Constraints bcs, DofType dof,
                      int numIndependentDofs, double time, std::string solver="EigenSparseLU");


} // namespace NuTo
