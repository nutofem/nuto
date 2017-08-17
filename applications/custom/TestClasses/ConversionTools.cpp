#include "ConversionTools.h"
#include "PrintTools.h"

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>



Epetra_CrsMatrix ConversionTools::convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, bool rAsGlobal, bool rFillComplete)
{
    int rowCount = rEigenMatrix.rows();
    int columnCount = rEigenMatrix.cols();

    if (rAsGlobal)
    {
        Epetra_Map rowMap(rowCount, 0, mComm);
        Epetra_Map columnMap(columnCount, 0, mComm);
        generateDefaultGlobalMaps(rowMap, columnMap, rowCount, columnCount);
        return convertEigen2EpetraCrsMatrix(rEigenMatrix, rowMap, columnMap, rFillComplete);
    }
    else
    {
        Epetra_Map rowMap(-1, rowCount, 0, mComm);
        Epetra_Map columnMap(-1, columnCount, 0, mComm);
        return convertEigen2EpetraCrsMatrix(rEigenMatrix, rowMap, columnMap, rFillComplete);
    }

}



Epetra_CrsMatrix ConversionTools::convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_Map rRowMap, Epetra_Map rColMap, bool rFillComplete)
{
    int rowCount = rEigenMatrix.rows();
    int columnCount = rEigenMatrix.cols();

    int numMyElements_Row = rRowMap.NumMyElements();
    int* myGlobalIndices_Row = rRowMap.MyGlobalElements();
    int numMyElements_Col = rColMap.NumMyElements();
    int* myGlobalIndices_Col = rColMap.MyGlobalElements();



    Epetra_CrsMatrix convMatrix(Copy, rRowMap, rColMap, 0);


    std::vector<int> localIndicesVector;
    std::vector<int> globalIndicesVector;
    int* localIndices;
    int* globalIndices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double value = 0.;

    int epetraErr = -1;

    for (int i = 0; i < numMyElements_Row; ++i)
//    for (int i = 0; i < rowCount; ++i)
    {
        numEntries = 0;
        for (int j = 0; j < numMyElements_Col; ++j)
//        for (int j = 0; j < columnCount; ++j)
        {
            value = rEigenMatrix.coeff(i, j);
            if (value != 0)
            {
                ++numEntries;
                entriesVector.push_back(value);
//                localIndicesVector.push_back(j);
                globalIndicesVector.push_back(myGlobalIndices_Col[j]);
            }
        }

        entries = &entriesVector[0];
//        localIndices = &localIndicesVector[0];
        globalIndices = &globalIndicesVector[0];

        convMatrix.InsertGlobalValues(myGlobalIndices_Row[i], numEntries, entries, globalIndices);
//        convMatrix.InsertMyValues(i, numEntries, entries, indices);

        entriesVector.clear();
        localIndicesVector.clear();
        globalIndicesVector.clear();
        numEntries = 0;
    }


    if (rFillComplete)
        convMatrix.FillComplete();

    return convMatrix;
}

//only for square matrices
Epetra_CrsMatrix ConversionTools::convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_CrsGraph rGraph, bool rFillComplete)
{
    Epetra_CrsMatrix convertedMatrix(Copy, rGraph);
    convertedMatrix.FillComplete();
    convertedMatrix.PutScalar(0.0);

    int rowCount = rEigenMatrix.rows();
    int columnCount = rEigenMatrix.cols();
    int* myGlobalIndices_Row = rGraph.Map().MyGlobalElements();
    int* myGlobalIndices_Col = rGraph.Map().MyGlobalElements();
    std::vector<int> globalIndicesVector;
    int* globalIndices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double value = 0.;


//    for (int i = 0; i < rowCount; ++i)
//    {
//        numEntries = 0;
//        for (int j = 0; j < columnCount; ++j)
//        {
//            value = rEigenMatrix.coeff(i, j);
//            if (value != 0)
//            {
//                ++numEntries;
//                entriesVector.push_back(value);
//                globalIndicesVector.push_back(myGlobalIndices_Col[j]);
//            }
//        }

//        entries = &entriesVector[0];
//        globalIndices = &globalIndicesVector[0];

//        convertedMatrix.SumIntoGlobalValues(myGlobalIndices_Row[i], numEntries, entries, globalIndices);

//        entriesVector.clear();
//        globalIndicesVector.clear();
//        numEntries = 0;
//    }

    for (int i = 0 ; i < rEigenMatrix.outerSize(); ++i)
    {
        numEntries = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(rEigenMatrix, i); it; ++it)
        {
            value = it.value();
            ++numEntries;
            entriesVector.push_back(value);
            globalIndicesVector.push_back(myGlobalIndices_Col[it.index()]);
        }

        entries = &entriesVector[0];
        globalIndices = &globalIndicesVector[0];

        convertedMatrix.SumIntoGlobalValues(myGlobalIndices_Row[i], numEntries, entries, globalIndices);

        entriesVector.clear();
        globalIndicesVector.clear();
        numEntries = 0;
    }

    if (rFillComplete)
        convertedMatrix.FillComplete();

    return convertedMatrix;
}


Epetra_Vector ConversionTools::convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, bool rAsGlobal)
{
    int length = rEigenVector.rows();

    if (rAsGlobal)
    {
        Epetra_Map map(length, 0, mComm);
        generateDefaultGlobalMap(map, length);

        return convertEigen2EpetraVector(rEigenVector, map);
    }
    else
    {
        Epetra_Map map(-1, length, 0, mComm);

        return convertEigen2EpetraVector(rEigenVector, map);
    }
}


Epetra_Vector ConversionTools::convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Epetra_Map rMap)
{
    int vectorSize = int(rEigenVector.rows());
    Epetra_Vector convertedVector(rMap);
    convertedVector.PutScalar(0.0);
    int* myGlobalIndices = rMap.MyGlobalElements();


    double* values = new double[1];
    int* globalIndices = new int[1];
    for (int i = 0; i < vectorSize; ++i)
    {
        values[0] = rEigenVector(i, 0);
        globalIndices[0] = myGlobalIndices[i];

//        convertedVector.SumIntoGlobalValues(1, values, globalIndices);
        convertedVector.ReplaceGlobalValues(1, values, globalIndices);
    }

    delete[] values;
    delete[] globalIndices;
//    delete myGlobalIndices;

    return convertedVector;

}


int* ConversionTools::map2Array_Int(std::map<int, int> rMap)
{
    int* newArray = new int[rMap.size()];

    for (int i = 0; i < rMap.size(); ++i)
    {
        newArray[i] = rMap[i];
    }

    return newArray;
}


std::map<int, int> ConversionTools::invertMap_int(std::map<int, int> rMap)
{
    std::map<int, int> newMap;

    for (int i = 0; i < rMap.size(); ++i)
    {
        newMap[rMap[i]] = i;
    }

    return newMap;
}

int ConversionTools::defineNumLocalElements(int rNumGlobalElements)
{
    int numProc = mComm.NumProc();
    int rank = mComm.MyPID();
    int numLocalElements = 0;

    double average = (double)rNumGlobalElements / numProc;
    if (rank < numProc-1)
        numLocalElements = std::floor(average);
    else
        numLocalElements = rNumGlobalElements - (std::floor(average)*(numProc-1));

    return numLocalElements;
}

void ConversionTools::generateDefaultGlobalMaps(Epetra_Map& rRowMap, Epetra_Map& rColumnMap, int rNumGlobalRows, int rNumGlobalColumns)
{
    int numLocalRows = defineNumLocalElements(rNumGlobalRows);

    Epetra_Map newRowMap(rNumGlobalRows, numLocalRows, 0, mComm);
    Epetra_Map newColumnMap(rNumGlobalColumns, rNumGlobalColumns, 0, mComm);

    rRowMap = newRowMap;
    rColumnMap = newColumnMap;
}

void ConversionTools::generateDefaultGlobalMap(Epetra_Map& rMap, int rNumGlobalElements)
{
    int numLocalElements = defineNumLocalElements(rNumGlobalElements);

    Epetra_Map newMap(rNumGlobalElements, numLocalElements, 0, mComm);

    rMap = newMap;
}




















