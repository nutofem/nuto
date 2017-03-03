#include "SolveSystem.h"
#include "mechanics/dofSubMatrixSolvers/SolverEigen.h"


BOOST_AUTO_TEST_CASE(SolverEigenSparseLU)
{
    NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>> s;
    SolveAndCheckSystem(s);
}

BOOST_AUTO_TEST_CASE(SolverEigenSparseQR)
{
    NuTo::SolverEigen<Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>> s;
    SolveAndCheckSystem(s);
}

BOOST_AUTO_TEST_CASE(SolverEigenSimplicialLDLT)
{
    NuTo::SolverEigen<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> s;
    SolveAndCheckSystem(s);
}
