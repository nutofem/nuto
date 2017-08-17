#include "../Benchmark.h"
#include "EpetraLinearProblemBenchmark.h"

#include <mpi.h>
#include <Epetra_MpiComm.h>

#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

#include "base/Timer.h"


const std::vector<int> numElements = {10, 10, 10};
const int numOfIterations = 1;


//BENCHMARK(Solver_Aztec_GlobalStructure, CG, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm);

//    while(runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_AztecOO("CG", true);
//    }
//}


BENCHMARK(Solver_Aztec_LocalStructure, CG, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

    while(runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveIterative_AztecOO("CG", true);
    }
}


//BENCHMARK(Solver_Aztec_GlobalStructure, CGS, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm);

//    while(runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_AztecOO("CGS", true);
//    }
//}


BENCHMARK(Solver_Aztec_LocalStructure, CGS, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

    while(runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveIterative_AztecOO("CGS", true);
    }
}


//BENCHMARK(Solver_Aztec_GlobalStructure, GMRES, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm);

//    while(runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_AztecOO("GMRES", true);
//    }
//}


BENCHMARK(Solver_Aztec_LocalStructure, GMRES, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

    while(runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveIterative_AztecOO("GMRES", true);
    }
}


//BENCHMARK(Solver_Aztec_GlobalStructure, BiCGStab, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm);

//    while(runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_AztecOO("BiCGStab", true);
//    }
//}


BENCHMARK(Solver_Aztec_LocalStructure, BiCGStab, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

    while(runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveIterative_AztecOO("BiCGStab", true);
    }
}


BENCHMARK(Solver_Belos_GlobalStructure, GMRES, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, true);

    while (runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveIterative_Belos("GMRES");
    }
}

BENCHMARK(Solver_Belos_LocalStructure, GMRES, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

    while (runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveIterative_Belos("GMRES");
    }
}

//BENCHMARK(Solver_Belos_GlobalStructure, GMRES_Precond, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, true);

//    while (runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_Belos("GMRES", true);
//    }
//}

//BENCHMARK(Solver_Belos_LocalStructure, GMRES_Precond, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

//    while (runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_Belos("GMRES", true);
//    }
//}


//BENCHMARK(Solver_Belos_GlobalStructure, CG, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, true);

//    while (runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_Belos("CG");
//    }
//}

BENCHMARK(Solver_Belos_LocalStructure, CG, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

    while (runner.KeepRunningIterations(numOfIterations))
    {
        epProblem.solveIterative_Belos("CG");
    }
}

//BENCHMARK(Solver_Belos_GlobalStructure, CG_Precond, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, true);

//    while (runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_Belos("CG", true);
//    }
//}

//BENCHMARK(Solver_Belos_LocalStructure, CG_Precond, runner)
//{
//    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements);
//    NuTo::Benchmark::EpetraLinearProblemBenchmark epProblem(s, comm, false);

//    while (runner.KeepRunningIterations(numOfIterations))
//    {
//        epProblem.solveIterative_Belos("CG", true);
//    }
//}



