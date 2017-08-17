#include "../Benchmark.h"
#include "EpetraLinearProblemBenchmark.h"

#include <mpi.h>
#include <Epetra_MpiComm.h>

#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

#include "base/Timer.h"


const std::vector<int> numElements = {5, 5, 5};
const int numOfIterations = 5;

BENCHMARK(Solver, Klu, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    int numProc = comm.NumProc();

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements); //, numProc);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm);

    while(runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveDirect("Klu");
    }
}


BENCHMARK(Solver, MUMPS, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm);


    while(runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveDirect("Mumps");
    }
}


BENCHMARK(Solver, Lapack, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm);


    while(runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveDirect("Lapack");
    }
}











