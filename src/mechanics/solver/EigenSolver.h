#pragma once
#include <string>
#include "Eigen/Sparse"

namespace NuTo
{

//! Solve a sparse linear system \f$A x = b\f$ made of Eigen types.
//! @param A Sparse matrix.
//! @param b Right hand side vector.
//! @param solver String representation of solver.
//!
//! Choose among the following solvers:
//!
//! Built-in direct solvers:
//! - `EigenSparseLU`
//! - `EigenSparseQR`
//! - `EigenSimplicialLLT`
//! - `EigenSimplicialLDLT`
//! 
//! Built-in iteratitive solvers:
//! - `EigenConjugateGradient`
//! - `EigenLeastSquaresConjugateGradient`
//! - `EigenBiCGSTAB.
//! 
//! and [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) solvers (need to be installed seperately):
//! - `EigenUmfPackLU`
//! - `EigenCholmodSupernodalLLT`
Eigen::VectorXd EigenSolver(Eigen::SparseMatrix<double> A, Eigen::VectorXd b, std::string solver);

} // namespace NuTo
