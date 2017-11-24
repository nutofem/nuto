#include <benchmark/benchmark.h>

#include "LinearElasticBenchmarkStructure.h"
#include "math/SparseMatrixCSR.h"

#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "mechanics/dofSubMatrixSolvers/SolverPardiso.h"
#include "mechanics/dofSubMatrixSolvers/SolverEigen.h"

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

// Conversion Benchmarks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
static void Convert_NuToToEigen(benchmark::State& state)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    for (auto _ : state)
        t.NuToToEigen();
}
BENCHMARK(Convert_NuToToEigen);


static void Convert_NuToToCSR(benchmark::State& state)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    for (auto _ : state)
        t.NuToToCSR();
}
BENCHMARK(Convert_NuToToCSR);


// Solver Benchmarks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <typename TSolver, typename... Params>
static void SolverBenchmark(benchmark::State& state, Params&&... params)
{
    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    TestProblem t(s);
    TSolver solver(std::forward<Params>(params)...);
    for (auto _ : state)
        t.Solve(solver);
}

static void Solver_MUMPS(benchmark::State& state)
{
    SolverBenchmark<NuTo::SolverMUMPS>(state, false);
}
BENCHMARK(Solver_MUMPS);

#ifdef HAVE_PARDISO
static void Solver_Pardiso(Solver_Pardiso)
{
    SolverBenchmark<NuTo::SolverPardiso>(state, 1, false);
}
BENCHMARK(Solver_Pardiso);
#endif /* HAVE_PARDISO */


static void Solver_EigenLU(benchmark::State& state)
{
    SolverBenchmark<NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>>(state);
}
BENCHMARK(Solver_EigenLU);

static void Solver_EigenLDLT(benchmark::State& state)
{
    SolverBenchmark<NuTo::SolverEigen<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>>(state);
}
BENCHMARK(Solver_EigenLDLT);


BENCHMARK_MAIN()
