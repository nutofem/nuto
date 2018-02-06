#include "MumpsSolver.h"
#include "base/Exception.h"
#include "MUMPSSupport"

using namespace NuTo;

template <typename TSolver>
Eigen::VectorXd MumpsSolve(Eigen::SparseMatrix<double> A, Eigen::VectorXd b)
{
    TSolver solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    return solver.solve(b);
}

Eigen::VectorXd NuTo::MumpsSolver(Eigen::SparseMatrix<double> A, Eigen::VectorXd b, std::string solver)
{
    if (solver == "MumpsLU")
        return MumpsSolve<Eigen::MUMPSLU<Eigen::SparseMatrix<double>>>(A, b);
    if (solver == "MumpsLDLT")
        return MumpsSolve<Eigen::MUMPSLDLT<Eigen::SparseMatrix<double>, Eigen::Upper>>(A, b);
    throw Exception("Unknown solver. Are you sure you spelled it correctly?");
}

