#include "Benchmark.h"
#include "LinearElasticBenchmarkStructure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"

#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "mechanics/dofSubMatrixSolvers/SolverEigen.h"
#include "mechanics/dofSubMatrixSolvers/SolverPardiso.h"

class TestProblem
{
public:
    TestProblem(NuTo::Benchmark::LinearElasticBenchmarkStructure& r)
        : matrix(r.GetStructure().GetDofStatus())
        , rhs(r.GetStructure().GetDofStatus())
    {
        r.SetupBCs();
        matrix = r.GetStructure().BuildGlobalHessian0().JJ;
        rhs = r.GetStructure().BuildGlobalInternalGradient().J;
    }

    void Solve(NuTo::SolverBase& rSolver)
    {
        rSolver.Solve(matrix, rhs);
    }

    void NuToToEigen()
    {
        matrix.ExportToEigenSparseMatrix();
    }

    void NuToToCSR()
    {
        matrix.ExportToCSR();
    }

private:
    NuTo::BlockSparseMatrix matrix;
    NuTo::BlockFullVector<double> rhs;
};

const std::vector<int> numElements = {10, 10, 10};

BENCHMARK(Convert, NuToToEigen, runner)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    while (runner.KeepRunningIterations(10))
    {
        t.NuToToEigen();
    }
}

BENCHMARK(Convert, NuToToCSR, runner)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    while (runner.KeepRunningIterations(10))
    {
        t.NuToToCSR();
    }
}

BENCHMARK(Solver, MUMPS, runner)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    while (runner.KeepRunningIterations(1))
    {
        NuTo::SolverMUMPS solver(false);
        t.Solve(solver);
    }
}

#ifdef HAVE_PARDISO
BENCHMARK(Solver, Pardiso, runner)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    while (runner.KeepRunningIterations(1))
    {
        NuTo::SolverPardiso solver(1, false);
        t.Solve(solver);
    }
}
#endif /* HAVE_PARDISO */

BENCHMARK(Solver, EigenLU, runner)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    while (runner.KeepRunningIterations(1))
    {
        NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>> solver;
        t.Solve(solver);
    }
}


BENCHMARK(Solver, EigenLDLT, runner)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    while (runner.KeepRunningIterations(1))
    {
        NuTo::SolverEigen<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>> solver;
        t.Solve(solver);
    }
}
