//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
#include "BoostUnitTest.h"

#include <mpi.h>

#include <Epetra_DataAccess.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack.h>
#include <Ifpack_AdditiveSchwarz.h>
//#include <Ifpack2_AdditiveSchwarz.hpp>
//#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include "Utilities.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


void compareVectors(Epetra_MultiVector rMVec1, Epetra_MultiVector rMVec2, int rMVecIndex1 = 0, int rMVecIndex2 = 0, double rTol = 1e-5)
{
    Epetra_Vector* vec1 = rMVec1(rMVecIndex1);
    Epetra_Vector* vec2 = rMVec2(rMVecIndex2);

    int n = vec1->MyLength();
    double* vals1 = vec1->Values();
    double* vals2 = vec2->Values();
    for (int i = 0; i < n; ++i)
    {
        BOOST_CHECK_CLOSE(vals1[i], vals2[i], rTol);
    }
}


#ifdef HAVE_MPI
void checkSolver(Epetra_MpiComm rComm)
#else
void checkSolver(Epetra_SerialComm rComm)
#endif
{
    int procID = rComm.MyPID();
    int errCode = 0;
    int numNonZerosPerRow = 3;
    Epetra_Map map(-1, numLocalElements, 0, rComm);
    Epetra_CrsMatrix A(Epetra_DataAccess::Copy, map, numNonZerosPerRow, false);
    Epetra_Vector rhs(map);
//    errCode = createLinearVector(rhs);
    errCode = rhs.ReplaceGlobalValue(0, 0, 1);
//    BOOST_CHECK_EQUAL(errCode, 0);
    Epetra_Vector lhs(map);
    Epetra_Vector checkVec(map);
    checkVec.PutScalar(1.);

    errCode = createMatrix(A);
    BOOST_CHECK_EQUAL(errCode, 0);
    errCode = A.FillComplete();
    BOOST_CHECK_EQUAL(errCode, 0);

    Epetra_LinearProblem prob(&A, &lhs, &rhs);

    if (procID == 0)
    {
        std::cout << "[Test_IterativeSolver]: PROBLEM DEFINITION" << std::endl;
        std::cout << "[Test_IterativeSolver]: System matrix" << std::endl;
    }
    A.Print(std::cout);

    if (procID == 0)
    {
        std::cout << "\n\n[Test_IterativeSolver]: Right hand side" << std::endl;
    }
    rhs.Print(std::cout);

    AztecOO az_solver(prob);
    az_solver.SetAztecOption(AZ_solver, AZ_gmres);
    az_solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    az_solver.SetAztecOption(AZ_kspace, 50);
    az_solver.SetAztecOption(AZ_output, AZ_all);
    az_solver.Iterate(1000, 1e-8);

    Epetra_MultiVector* solution;
    solution = prob.GetLHS();
    compareVectors(*solution, checkVec);

    if (procID == 0)
        std::cout << "\n\n[Test_IterativeSolver]: Solution" << std::endl;
    solution->Print(std::cout);
}



BOOST_AUTO_TEST_CASE(IterativeSolver)
{
    MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    checkSolver(comm);

    MPI_Finalize();
}
































