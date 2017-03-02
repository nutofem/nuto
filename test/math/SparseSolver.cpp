//
// Created by phuschke on 3/1/17.
//

#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Sparse>
#include "BoostUnitTest.h"
#include "mechanics/dofSubMatrixSolvers/SolverEigen.h"
#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSR.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/dofSubMatrixStorage/BlockStorageBase.h"

using NuTo::Node::eDof;

void SolveAndCheckSystem(NuTo::SolverBase& solver,
                         const NuTo::BlockSparseMatrix& matrix,
                         const NuTo::BlockFullVector<double>& rhs,
                         const NuTo::BlockFullVector<double>& expectedSolution)
{
    BOOST_CHECK( (solver.Solve(matrix, rhs).Export() - expectedSolution.Export()).isMuchSmallerThan(1.e-6,1.e-1) );
}

void TestDofSubMatrixSolvers(const NuTo::BlockSparseMatrix& matrix,
                             const NuTo::BlockFullVector<double>& rhs,
                             const NuTo::BlockFullVector<double>& expectedSolution)
{

    // direct solvers
    NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>> solverLU;
    SolveAndCheckSystem(solverLU, matrix, rhs, expectedSolution);

    NuTo::SolverEigen<Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>> solverQR;
    SolveAndCheckSystem(solverQR, matrix, rhs, expectedSolution);

    NuTo::SolverEigen<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> solverLDLT;
    SolveAndCheckSystem(solverLDLT, matrix, rhs, expectedSolution);

    NuTo::SolverMUMPS solverMUMPS;
    SolveAndCheckSystem(solverMUMPS, matrix, rhs, expectedSolution);

    /// \todo add pardiso, pardiso mkl, external eigen wrappers

}

BOOST_AUTO_TEST_CASE(dofSubMatrixSolvers)
{

    const int submatrixSize = 2;
    std::set<NuTo::Node::eDof> dofTypes;
    dofTypes.insert(eDof::DISPLACEMENTS);
    dofTypes.insert(eDof::CRACKPHASEFIELD);
    
    NuTo::DofStatus dofStatus;
    dofStatus.SetDofTypes(dofTypes);
    dofStatus.SetActiveDofTypes(dofTypes);
    dofStatus.SetNumActiveDofs(eDof::DISPLACEMENTS, submatrixSize);
    dofStatus.SetNumActiveDofs(eDof::CRACKPHASEFIELD, submatrixSize);

    NuTo::BlockFullVector<double> rhs(dofStatus);
    rhs[eDof::DISPLACEMENTS]    = Eigen::Vector2d::Constant(1.);
    rhs[eDof::CRACKPHASEFIELD]  = Eigen::Vector2d::Constant(1.);

    NuTo::BlockSparseMatrix matrix(dofStatus);

    matrix(eDof::DISPLACEMENTS,   eDof::DISPLACEMENTS).Resize(submatrixSize,submatrixSize);
    matrix(eDof::DISPLACEMENTS,   eDof::CRACKPHASEFIELD).Resize(submatrixSize,submatrixSize);
    matrix(eDof::CRACKPHASEFIELD, eDof::DISPLACEMENTS).Resize(submatrixSize,submatrixSize);
    matrix(eDof::CRACKPHASEFIELD, eDof::CRACKPHASEFIELD).Resize(submatrixSize,submatrixSize);

    matrix(eDof::DISPLACEMENTS, eDof::DISPLACEMENTS).AddValue(0,0,2.);
    matrix(eDof::DISPLACEMENTS, eDof::DISPLACEMENTS).AddValue(1,1,2.);

    matrix(eDof::CRACKPHASEFIELD, eDof::CRACKPHASEFIELD).AddValue(0,0,2.);
    matrix(eDof::CRACKPHASEFIELD, eDof::CRACKPHASEFIELD).AddValue(1,1,2.);

    NuTo::BlockFullVector<double> expectedSolution(dofStatus);
    expectedSolution[eDof::DISPLACEMENTS]    = Eigen::Vector2d::Constant(0.5);
    expectedSolution[eDof::CRACKPHASEFIELD]  = Eigen::Vector2d::Constant(0.5);

    TestDofSubMatrixSolvers(matrix, rhs, expectedSolution);

}