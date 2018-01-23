#include "EigenSolver.h"
#include "base/Exception.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

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
    if (solver == "EigenSparseLU")
        return EigenSparseSolve<Eigen::SparseLU<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "EigenSparseQR")
        return EigenSparseSolve<Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>(A, b);
    throw Exception("Que?");
}

} // namespace NuTo
