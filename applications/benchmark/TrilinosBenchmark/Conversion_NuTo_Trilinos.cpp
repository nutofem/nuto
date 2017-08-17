#include "../Benchmark.h"
#include "../LinearElasticBenchmarkStructure.h"

#include <mpi.h>
#include <Epetra_MpiComm.h>

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

#include "../../applications/custom/TestClasses/ConversionTools.h"


const std::vector<int> numElements = {10, 10, 10};


BENCHMARK(Conversion_NuTo_Trilinos, BlockSparse2EpetraCrs, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    int numProc = comm.NumProc();

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements); //, numProc);
    ConversionTools converter(comm);
    Eigen::SparseMatrix<double> mat_eigen = s.GetStructure().BuildGlobalHessian0().JJ.ExportToEigenSparseMatrix();
    Epetra_CrsMatrix* mat_epetra;

    while(runner.KeepRunningIterations(10))
    {
        mat_epetra = new Epetra_CrsMatrix(converter.convertEigen2EpetraCrsMatrix(mat_eigen));
    }
}



BENCHMARK(Conversion_NuTo_Trilinos, BlockSparse2EpetraCrs_global, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    int numProc = comm.NumProc();

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements); //, numProc);
    ConversionTools converter(comm);
    Eigen::SparseMatrix<double> mat_eigen = s.GetStructure().BuildGlobalHessian0().JJ.ExportToEigenSparseMatrix();
    Epetra_CrsMatrix* mat_epetra;

    while(runner.KeepRunningIterations(10))
    {
        mat_epetra = new Epetra_CrsMatrix(converter.convertEigen2EpetraCrsMatrix(mat_eigen, true));
    }
}


BENCHMARK(Conversion_NuTo_Trilinos, BlockVector2EpetraVector, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    int numProc = comm.NumProc();

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements); //, numProc);
    ConversionTools converter(comm);
    Eigen::VectorXd vec_eigen = s.GetStructure().BuildGlobalInternalGradient().J.Export();
    Epetra_Vector* vec_epetra;

    while(runner.KeepRunningIterations(10))
    {
        vec_epetra = new Epetra_Vector(converter.convertEigen2EpetraVector(vec_eigen));
    }
}


BENCHMARK(Conversion_NuTo_Trilinos, BlockVector2EpetraVector_global, runner)
{
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    int numProc = comm.NumProc();

    NuTo::Benchmark::LinearElasticBenchmarkStructure s(numElements); //, numProc);
    ConversionTools converter(comm);
    Eigen::VectorXd vec_eigen = s.GetStructure().BuildGlobalInternalGradient().J.Export();
    Epetra_Vector* vec_epetra;

    while(runner.KeepRunningIterations(10))
    {
        vec_epetra = new Epetra_Vector(converter.convertEigen2EpetraVector(vec_eigen, true));
    }
}
