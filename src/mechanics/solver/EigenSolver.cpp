#include "EigenSolver.h"
#include "base/Exception.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

#ifdef HAVE_SUITESPARSE
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#endif

namespace NuTo
{

template <typename TSolver>
Eigen::VectorXd EigenSparseSolve(Eigen::SparseMatrix<double> A, Eigen::VectorXd b)
{
    TSolver solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    return solver.solve(b);
}

Eigen::VectorXd EigenSolver(Eigen::SparseMatrix<double> A, Eigen::VectorXd b, std::string solver)
{
    // direct solvers
    if (solver == "EigenSparseLU")
        return EigenSparseSolve<Eigen::SparseLU<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenSparseQR")
        return EigenSparseSolve<Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>(A, b);
    if (solver == "EigenSimplicialLLT")
        return EigenSparseSolve<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenSimplicialLDLT")
        return EigenSparseSolve<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>(A, b);
    // iterative solvers
    if (solver == "EigenConjugateGradient")
        return EigenSparseSolve<Eigen::ConjugateGradient<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenLeastSquaresConjugateGradient")
        return EigenSparseSolve<Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenBiCGSTAB")
        return EigenSparseSolve<Eigen::BiCGSTAB<Eigen::SparseMatrix<double>>>(A, b);
    // external solvers
#ifdef HAVE_SUITESPARSE
    if (solver == "EigenUmfPackLU")
        return EigenSparseSolve<Eigen::UmfPackLU<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenCholmodSupernodalLLT")
        return EigenSparseSolve<Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>>(A, b);
#else
    if (solver == "EigenUmfPackLU" or solver == "EigenCholmodSupernodalLLT")
        throw Exception("NuTo has not been compiled with SuiteSparse.");
#endif
    throw Exception("Unknown solver. Are you sure you spelled it correctly?");
}

} // namespace NuTo
