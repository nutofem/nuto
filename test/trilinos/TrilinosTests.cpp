//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
#include "BoostUnitTest.h"

#include <mpi.h>

#include <Epetra_DataAccess.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Export.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Amesos_ConfigDefs.h>
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
#include "../../applications/custom/TestClasses/ConversionTools.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

constexpr int numLocalElements = 5;

int createGraph(int rSetNonZerosPerRow, Epetra_CrsGraph& rGraph)
{
    int procID = rGraph.Comm().MyPID();
//    Epetra_CrsGraph graph(Epetra_DataAccess::Copy, rMap, rNonZerosPerRow, rSetStaticProfile);

    int errCode = 0;
    int globRow = 0;
    int numIndices = 1;
    int* indices = new int[1];
    for (int i = 0; i < rSetNonZerosPerRow; ++i)
    {
        globRow = i + procID*numLocalElements;
        numIndices = i+1;
        indices = new int[numIndices];
        for (int j = 0; j < numIndices; ++j)
            indices[j] = j + procID*numLocalElements;

        errCode = rGraph.InsertGlobalIndices(globRow, numIndices, indices);
        if (errCode != 0)
            break;
    }

    delete[] indices;

    return errCode;
}


int createGraph_2(Epetra_CrsGraph& rGraph)
{
    int mySize = rGraph.Map().NumMyElements();
    int maxSize = rGraph.Map().MaxAllGID();
    int* myGlobInds = rGraph.Map().MyGlobalElements();
    int errCode = 0;
    int globRow = 0;
    int numIndices = 1;
    int* indices;
    for (int i = 0; i < mySize; ++i)
    {
        globRow = myGlobInds[i];
        if (globRow == 0)
        {
            numIndices = 2;
            indices = new int[numIndices];
            indices[0] = 0;
            indices[1] = 1;
        }
        else if (globRow == maxSize)
        {
            numIndices = 2;
            indices = new int[numIndices];
            indices[0] = maxSize-1;
            indices[1] = maxSize;
        }
        else
        {
            numIndices = 3;
            indices = new int[numIndices];
            indices[0] = globRow-1;
            indices[1] = globRow;
            indices[2] = globRow+1;
        }

        errCode = rGraph.InsertGlobalIndices(globRow, numIndices, indices);
        if (errCode != 0)
            break;
    }

    delete[] indices;
//    delete myGlobInds;

    return errCode;
}


int createMatrix(Epetra_CrsMatrix& rMatrix, bool graphProvided = false)
{
    int errCode = 0;
    int rowCount = rMatrix.NumMyRows();
    int* myGlobInd = rMatrix.Map().MyGlobalElements();
    int maxCount = rMatrix.Map().MaxAllGID();
    int* globIndices;
    int numEntries = 0;
    double* entries;

    for (int i = 0; i < rowCount; ++i)
    {
        if (myGlobInd[i] == 0)
        {
            numEntries = 2;
            entries = new double[numEntries];
            globIndices = new int[numEntries];
            entries[0] = 2;
            entries[1] = -1;
            globIndices[0] = 0;
            globIndices[1] = 1;
        }
        else if (myGlobInd[i] == maxCount)
        {
            numEntries = 2;
            entries = new double[numEntries];
            globIndices = new int[numEntries];
            entries[0] = -1;
            entries[1] = 1;
            globIndices[0] = maxCount-1;
            globIndices[1] = maxCount;
        }
        else
        {
            numEntries = 3;
            entries = new double[numEntries];
            globIndices = new int[numEntries];
            entries[0] = -1;
            entries[1] = 2;
            entries[2] = -1;
            globIndices[0] = myGlobInd[i]-1;
            globIndices[1] = myGlobInd[i];
            globIndices[2] = myGlobInd[i]+1;
        }

        if (graphProvided)
            errCode = rMatrix.SumIntoGlobalValues(myGlobInd[i], numEntries, entries, globIndices);
        else
            errCode = rMatrix.InsertGlobalValues(myGlobInd[i], numEntries, entries, globIndices);

        if (errCode != 0)
            break;
    }

//    delete myGlobInd;
    delete[] globIndices;
    delete[] entries;

    return errCode;
}


int createLinearVector(Epetra_Vector& rVector)
{
    int errCode = 0;
    int n = rVector.MyLength();
    int* globIndices = rVector.Map().MyGlobalElements();

    for (int i = 0; i < n; ++i)
    {
        errCode = rVector.ReplaceGlobalValue(globIndices[i], 0, globIndices[i]);
        if (errCode != 0)
            break;
    }

//    delete globIndices;

    return errCode;
}


#ifdef HAVE_MPI
void checkEpetraGraphCreation(Epetra_MpiComm rComm)
#else
void checkEpetraGraphCreation(Epetra_SerialComm rComm)
#endif
{
    int procCount = rComm.NumProc();
    int procID = rComm.MyPID();
//    int pcount  = 4;
//    int pid = 1;
//    int numLocalElements = 5;
    int numGlobalElements = procCount * numLocalElements;

    int* specificMapping = new int[numLocalElements];
    for (int i = 0; i < numLocalElements; ++i)
    {
        specificMapping[i] = procCount*i+procID;
    }

    Epetra_Map linearMap(-1, numLocalElements, 0, rComm);
    Epetra_Map specificMap(-1, numLocalElements, specificMapping, 0, rComm);

    int maxLID = numLocalElements-1;
    int maxGID = numGlobalElements-1;
    int myMaxGID = procCount*(numLocalElements-1) + procID;
//    BOOST_CHECK_EQUAL(maxGID, myMaxGID);
    BOOST_CHECK_EQUAL(specificMap.MaxLID(), maxLID);
    BOOST_CHECK_EQUAL(specificMap.MaxMyGID(), myMaxGID);
    BOOST_CHECK_EQUAL(specificMap.MaxAllGID(), maxGID);

    int nonZerosPerRow = 3;
    bool setStaticProfile = false;
    int errCode = 0;

    Epetra_CrsGraph linearGraph_staticFalse(Epetra_DataAccess::Copy, linearMap, nonZerosPerRow, setStaticProfile);
    errCode = createGraph(nonZerosPerRow, linearGraph_staticFalse);
    BOOST_CHECK_EQUAL(errCode, 0);
    linearGraph_staticFalse.FillComplete();
    BOOST_CHECK_CLOSE(double(linearGraph_staticFalse.NumMyNonzeros()), nonZerosPerRow*(nonZerosPerRow+1)/2, 0.001);

    Epetra_CrsGraph linearGraph_staticFalse_augment(Epetra_DataAccess::Copy, linearMap, nonZerosPerRow, setStaticProfile);
    errCode = createGraph(nonZerosPerRow+1, linearGraph_staticFalse_augment);
    BOOST_CHECK_EQUAL(errCode, 0);
    linearGraph_staticFalse_augment.FillComplete();
    BOOST_CHECK_CLOSE(double(linearGraph_staticFalse_augment.NumMyNonzeros()), (nonZerosPerRow+1)*(nonZerosPerRow+2)/2, 0.001);

    setStaticProfile = true;
    Epetra_CrsGraph linearGraph_staticTrue_augment(Epetra_DataAccess::Copy, linearMap, nonZerosPerRow, setStaticProfile);
    errCode = createGraph(nonZerosPerRow+1, linearGraph_staticTrue_augment);
    BOOST_CHECK_LT(errCode, 0); //if static profile is set to true, it has to occure some error
    linearGraph_staticTrue_augment.FillComplete();
    BOOST_CHECK_CLOSE(double(linearGraph_staticTrue_augment.NumMyNonzeros()), nonZerosPerRow*(nonZerosPerRow+1)/2, 0.001);
}


#ifdef HAVE_MPI
void checkEpetraExport(Epetra_MpiComm rComm)
#else
void checkEpetraExport(Epetra_SerialComm rComm)
#endif
{
    int procCount = rComm.NumProc();
    int procID = rComm.MyPID();
    int offset = 0;
    if (procID > 0)
    {
        offset = 1;
    }

    int* overlapMapping = new int[numLocalElements + offset];

    for (int i = 0; i < numLocalElements; ++i)
    {
        overlapMapping[i] = procID*numLocalElements + i;
    }

    for (int k = 0; k < offset; ++k)
    {
        overlapMapping[numLocalElements+k] = procID*numLocalElements-offset+k;
    }

    Epetra_Map specificMap(-1, numLocalElements, 0, rComm);
    Epetra_Map overlapMap(-1, numLocalElements + offset, overlapMapping, 0, rComm);
    Epetra_Export exporter(overlapMap, specificMap);

    int errCode = 0;
    int nonZerosPerRow = 3;
    bool setStaticProfile = false;
    Epetra_CrsGraph specificGraph(Epetra_DataAccess::Copy, specificMap, nonZerosPerRow, setStaticProfile);
    Epetra_CrsGraph overlapGraph(Epetra_DataAccess::Copy, overlapMap, nonZerosPerRow, setStaticProfile);

    errCode = createGraph_2(overlapGraph);

    errCode = overlapGraph.FillComplete();
    errCode = specificGraph.Export(overlapGraph, exporter, Epetra_CombineMode::Insert);
    errCode = specificGraph.FillComplete();

    Epetra_CrsMatrix globMatrix(Epetra_DataAccess::Copy, specificGraph);
    Epetra_CrsMatrix locMatrix(Epetra_DataAccess::Copy, overlapGraph);
    errCode = createMatrix(locMatrix, true);

    Epetra_Vector globVector(specificMap);
    globVector.PutScalar(0.);
    Epetra_Vector locVector(overlapMap);
    locVector.PutScalar(double(procID));

    Epetra_CombineMode* combineModes = new Epetra_CombineMode[3];
    combineModes[0] = Epetra_CombineMode::Insert;
    combineModes[1] = Epetra_CombineMode::Add;
    combineModes[2] = Epetra_CombineMode::Average;

    double checkValue = 2.;
    double* rowValues = new double[nonZerosPerRow];
    int* colIndices = new int[nonZerosPerRow];
    int numEntries = 0;
    //combine mode "Average" not supported for matrices
    for (int c = 0; c < 2; ++c)
    {
        errCode = globMatrix.Export(locMatrix, exporter, combineModes[c]);
        BOOST_CHECK_EQUAL(errCode, 0);
        errCode = globMatrix.ExtractGlobalRowCopy((procID+1)*numLocalElements-1, nonZerosPerRow, numEntries, rowValues, colIndices);
        BOOST_CHECK_EQUAL(errCode, 0);
        checkValue = (combineModes[c] == Epetra_CombineMode::Insert ? 2. : 4.);
        if (procID < procCount-1)
            BOOST_CHECK_EQUAL(rowValues[1], checkValue);
    }
    delete[] rowValues;
    delete[] colIndices;

    for (int c = 0; c < 3; ++c)
    {
        errCode = globVector.Export(locVector, exporter, combineModes[c]);
        BOOST_CHECK_EQUAL(errCode, 0);
        switch (combineModes[c]) {
        case Epetra_CombineMode::Insert:
            checkValue = procID+1;
            break;
        case Epetra_CombineMode::Add:
            checkValue = 2*procID+1;
            break;

        case Epetra_CombineMode::Average:
            checkValue = (2*procID+1)/2;
            break;
        default:
            checkValue = 0;
            break;
        }
    }

    delete[] combineModes;
}


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

    AztecOO az_solver(prob);
    az_solver.SetAztecOption(AZ_solver, AZ_gmres);
    az_solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    az_solver.SetAztecOption(AZ_kspace, 50);
    az_solver.SetAztecOption(AZ_output, AZ_none);
    az_solver.Iterate(1000, 1e-8);
    compareVectors(*prob.GetLHS(), checkVec);

    int numSolvers = 3;
    Amesos Factory;
    std::string* solverTypes = new std::string[numSolvers];
    solverTypes[0] = "Klu";
    solverTypes[1] = "Lapack";
    solverTypes[2] = "Mumps";

    bool solverAvailable = true;

    for (int i = 0; i < numSolvers; ++i)
    {
        solverAvailable = Factory.Query(solverTypes[i]);
        BOOST_CHECK_EQUAL(solverAvailable, true);
    }

    Teuchos::ParameterList params;
    params.set("PrintStatus", true);
    params.set("PrintTiming", true);
    params.set("MaxProcs", -3); //all processes in communicator will be used

    Amesos_BaseSolver* am_solver;

    for (int i = 0; i < numSolvers; ++i)
    {
        am_solver = Factory.Create(solverTypes[i], prob);
        am_solver->SetParameters(params);
        am_solver->Solve();

        compareVectors(*prob.GetLHS(), checkVec);
    }

}



BOOST_AUTO_TEST_CASE(EpetraUtilities)
{
    MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    checkEpetraGraphCreation(comm);

//    checkEpetraExport(comm);

    checkSolver(comm);

    MPI_Finalize();
}



BOOST_AUTO_TEST_CASE(ConversionUtilities)
{
//    MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
//             &boost::unit_test::framework::master_test_suite().argv);

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    checkEpetraExport(comm);


//    MPI_Finalize();
}
































