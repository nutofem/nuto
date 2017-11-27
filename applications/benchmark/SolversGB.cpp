#include <benchmark/benchmark.h>

#include "LinearElasticBenchmarkStructure.h"
#include "math/SparseMatrixCSR.h"

#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "mechanics/dofSubMatrixSolvers/SolverPardiso.h"
#include "mechanics/dofSubMatrixSolvers/SolverEigen.h"


// Setup Test Structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

// Benchmark Fixture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class fixture : public benchmark::Fixture
{
    const std::vector<int> numElements = {10, 10, 10};
    NuTo::Benchmark::LinearElasticBenchmarkStructure s = {numElements};

public:
    TestProblem t = {s};

    template <typename TSolver = void, typename... Params>
    void SolverBenchmark(benchmark::State& state, Params&&... params)
    {
        TSolver solver = {std::forward<Params>(params)...};
        for (auto _ : state)
            t.Solve(solver);
    }
};


// Conversion Benchmarks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BENCHMARK_F(fixture, Convert_NuToToEigen)(benchmark::State& state)
{
    for (auto _ : state)
        t.NuToToEigen();
}


BENCHMARK_F(fixture, Convert_NuToToCSR)(benchmark::State& state)
{
    for (auto _ : state)
        t.NuToToCSR();
}


// Solver Benchmarks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BENCHMARK_F(fixture, Solver_MUMPS)(benchmark::State& state)
{
    SolverBenchmark<NuTo::SolverMUMPS>(state, false);
}


#ifdef HAVE_PARDISO
BENCHMARK_F(fixture, Solver_Pardiso)(Solver_Pardiso)
{
    SolverBenchmark<NuTo::SolverPardiso>(state, 1, false);
}
#endif /* HAVE_PARDISO */


BENCHMARK_F(fixture, Solver_EigenLU)(benchmark::State& state)
{
    SolverBenchmark<NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>>(state);
}

BENCHMARK_F(fixture, Solver_EigenLDLT)(benchmark::State& state)
{
    SolverBenchmark<NuTo::SolverEigen<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>>(state);
}


BENCHMARK_MAIN()
