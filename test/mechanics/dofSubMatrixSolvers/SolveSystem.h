#include "BoostUnitTest.h"
#include "mechanics/dofSubMatrixSolvers/SolverBase.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/nodes/NodeEnum.h"


struct TestProblem
{
    TestProblem()
        : rhs(dofStatus)
        , matrix(dofStatus)
        , expectedSolution(dofStatus)
    {
        using NuTo::Node::eDof;
        const int submatrixSize = 2;
        std::set<NuTo::Node::eDof> dofTypes({eDof::DISPLACEMENTS, eDof::CRACKPHASEFIELD});

        dofStatus.SetDofTypes(dofTypes);
        dofStatus.SetActiveDofTypes(dofTypes);
        dofStatus.SetNumActiveDofs(eDof::DISPLACEMENTS, submatrixSize);
        dofStatus.SetNumActiveDofs(eDof::CRACKPHASEFIELD, submatrixSize);

        rhs.AllocateSubvectors();
        rhs[eDof::DISPLACEMENTS] = Eigen::Vector2d::Constant(1.);
        rhs[eDof::CRACKPHASEFIELD] = Eigen::Vector2d::Constant(1.);

        matrix.AllocateSubmatrices();
        matrix(eDof::DISPLACEMENTS, eDof::DISPLACEMENTS).Resize(submatrixSize, submatrixSize);
        matrix(eDof::DISPLACEMENTS, eDof::CRACKPHASEFIELD).Resize(submatrixSize, submatrixSize);
        matrix(eDof::CRACKPHASEFIELD, eDof::DISPLACEMENTS).Resize(submatrixSize, submatrixSize);
        matrix(eDof::CRACKPHASEFIELD, eDof::CRACKPHASEFIELD).Resize(submatrixSize, submatrixSize);

        matrix(eDof::DISPLACEMENTS, eDof::DISPLACEMENTS).AddValue(0, 0, 2.);
        matrix(eDof::DISPLACEMENTS, eDof::DISPLACEMENTS).AddValue(1, 1, 2.);

        matrix(eDof::CRACKPHASEFIELD, eDof::CRACKPHASEFIELD).AddValue(0, 0, 2.);
        matrix(eDof::CRACKPHASEFIELD, eDof::CRACKPHASEFIELD).AddValue(1, 1, 2.);

        expectedSolution.AllocateSubvectors();
        expectedSolution[eDof::DISPLACEMENTS] = Eigen::Vector2d::Constant(0.5);
        expectedSolution[eDof::CRACKPHASEFIELD] = Eigen::Vector2d::Constant(0.5);
    }


    NuTo::DofStatus dofStatus;
    NuTo::BlockFullVector<double> rhs;
    NuTo::BlockSparseMatrix matrix;
    NuTo::BlockFullVector<double> expectedSolution;
};

void SolveAndCheckSystem(NuTo::SolverBase& solver)
{
    TestProblem p;
    BOOST_CHECK((solver.Solve(p.matrix, p.rhs).Export() - p.expectedSolution.Export()).isMuchSmallerThan(1.e-6, 1.e-1));
}
