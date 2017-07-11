//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
#include "BoostUnitTest.h"

#include <mpi.h>

#include <Epetra_DataAccess.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Amesos_ConfigDefs.h>

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
        std::cout << "[Test_DirectSolver]: PROBLEM DEFINITION" << std::endl;
        std::cout << "[Test_DirectSolver]: System matrix" << std::endl;
    }
    A.Print(std::cout);

    if (procID == 0)
    {
        std::cout << "\n\n[Test_DirectSolver]: Right hand side" << std::endl;
    }
    rhs.Print(std::cout);

    int numSolvers = 3;
    Amesos Factory;
    std::string* solverTypes = new std::string[numSolvers];
    solverTypes[0] = "Klu";
    solverTypes[1] = "Lapack";
    solverTypes[2] = "Mumps";
    bool* solverAvail = new bool[numSolvers];

    bool solverAvailable = true;

    for (int i = 0; i < numSolvers; ++i)
    {
        solverAvailable = Factory.Query(solverTypes[i]);
        solverAvail[i] = solverAvailable;
//        BOOST_CHECK_EQUAL(solverAvailable, true);
    }

    Teuchos::ParameterList params;
    params.set("PrintStatus", true);
    params.set("PrintTiming", true);
    params.set("MaxProcs", -3); //all processes in communicator will be used

    Amesos_BaseSolver* am_solver;
    Epetra_MultiVector* solution;
    for (int i = 0; i < numSolvers; ++i)
    {
        if (solverAvail[i])
        {
            am_solver = Factory.Create(solverTypes[i], prob);
            am_solver->SetParameters(params);
            am_solver->Solve();

            solution = prob.GetLHS();
            compareVectors(*solution, checkVec);

            if (procID == 0)
                std::cout << "\n\n[Test_DirectSolver]: Solution with '" << solverTypes[i] << "':" << std::endl;
            solution->Print(std::cout);
        }
        else
        {
            if (procID == 0)
                std::cout << "\n\n[Test_DirectSolver]: Solver '" << solverTypes[i] << "' not available!" << std::endl;
        }
    }

}



BOOST_AUTO_TEST_CASE(DirectSolver)
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

































