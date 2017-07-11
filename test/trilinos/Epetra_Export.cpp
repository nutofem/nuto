//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
#include "BoostUnitTest.h"

#include <mpi.h>

#include <Epetra_DataAccess.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Export.h>

#include "Utilities.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


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
    std::string* combineModesString = new std::string[3];
    combineModes[0] = Epetra_CombineMode::Insert;
    combineModes[1] = Epetra_CombineMode::Add;
    combineModes[2] = Epetra_CombineMode::Average;
    combineModesString[0] = "Insert";
    combineModesString[1] = "Add";
    combineModesString[2] = "Average";

    if (procID == 0)
    {
        std::cout << "[Test_EpetraExport]: Local matrices" << std::endl;
    }
    locMatrix.Print(std::cout);

    double checkValue = 2.;
    double* rowValues = new double[nonZerosPerRow];
    int* colIndices = new int[nonZerosPerRow];
    int numEntries = 0;
    //combine mode "Average" not supported for matrices
    for (int c = 0; c < 2; ++c)
    {
        errCode = globMatrix.Export(locMatrix, exporter, combineModes[c]);
        BOOST_CHECK_EQUAL(errCode, 0);
        if (procID == 0)
        {
            std::cout << "\n\n[Test_EpetraExport]: Global matrix with combine mode '" << combineModesString[c] << "'" << std::endl;
        }
        globMatrix.Print(std::cout);
        errCode = globMatrix.ExtractGlobalRowCopy((procID+1)*numLocalElements-1, nonZerosPerRow, numEntries, rowValues, colIndices);
        BOOST_CHECK_EQUAL(errCode, 0);
        checkValue = (combineModes[c] == Epetra_CombineMode::Insert ? 2. : 4.);
        if (procID < procCount-1)
            BOOST_CHECK_EQUAL(rowValues[1], checkValue);
    }
    delete[] rowValues;
    delete[] colIndices;

    if (procID == 0)
    {
        std::cout << "\n\n[Test_EpetraExport]: Local vectors" << std::endl;
    }
    locVector.Print(std::cout);

    for (int c = 0; c < 3; ++c)
    {
        errCode = globVector.Export(locVector, exporter, combineModes[c]);
        BOOST_CHECK_EQUAL(errCode, 0);
        if (procID == 0)
        {
            std::cout << "\n\n[Test_EpetraExport]: Global vector with combine mode '" << combineModesString[c] << "'" << std::endl;
        }
        globVector.Print(std::cout);
        switch (combineModes[c]) {
        case Epetra_CombineMode::Insert:
            checkValue = procID+1;
            break;
        case Epetra_CombineMode::Add:
            checkValue = 2*procID+1;
            break;

        case Epetra_CombineMode::Average:
            checkValue = (2*procID+1)/2.;
            break;
        default:
            checkValue = 0;
            break;
        }

        if (procID < procCount-1)
            BOOST_CHECK_EQUAL(globVector[numLocalElements-1], checkValue);
    }

    delete[] combineModes;
}



BOOST_AUTO_TEST_CASE(EpetraExport)
{
    MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    checkEpetraExport(comm);

    MPI_Finalize();
}
































