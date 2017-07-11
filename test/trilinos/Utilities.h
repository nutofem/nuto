//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
//#include "BoostUnitTest.h"

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


int createInverseLinearMap(Epetra_Map& rMap)
{
    int* mapping = new int[numLocalElements];
    int procID = rMap.Comm().MyPID();
    int maxValue = (procID+1)*numLocalElements-1;

    for (int i = 0; i < numLocalElements; ++i)
    {
        mapping[i] = maxValue - i;
    }

    Epetra_Map inverseMap(-1, numLocalElements, mapping, 0, rMap.Comm());

    rMap = inverseMap;

    delete[] mapping;
    return 0;
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
































