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
//! ## Built-in
//! ### Direct solvers:
//! - `EigenSparseLU`
//! - `EigenSparseQR`
//! - `EigenSimplicialLLT`
//! - `EigenSimplicialLDLT`
//!
//! ### Iteratitive solvers:
//! - `EigenConjugateGradient`
//! - `EigenLeastSquaresConjugateGradient`
//! - `EigenBiCGSTAB`
//!
//! ## External (need to be installed seperately)
//! ### [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) solvers:
//! - `SuiteSparseLU`
//! - `SuiteSparseSupernodalLLT`
//!
//! ### [MUMPS](http://mumps.enseeiht.fr/) solvers:
//! - `MumpsLU`
//! - `MumpsLDLT`
Eigen::VectorXd EigenSparseSolve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, std::string solver);

//! Solver usable by NewtonRaphson::Solve(...)
class EigenSparseSolver
{
public:
    EigenSparseSolver(std::string solver);
    Eigen::VectorXd Solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) const;

private:
    std::string mSolver;
};

} // namespace NuTo
