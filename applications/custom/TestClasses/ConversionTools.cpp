#include "ConversionTools.h"

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <eigen3/Eigen/Core>

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



Epetra_CrsMatrix ConversionTools::convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_Map rRowMap, Epetra_Map rColMap)
{
    int columnCount = rEigenMatrix.cols();

    int numMyElements_Row = rRowMap.NumMyElements();
//    int* myGlobalNodeIndices_Row = rRowMap.MyGlobalElements();


    Epetra_CrsMatrix convMatrix(Copy, rRowMap, rColMap, 0);


    std::vector<int> indicesVector;
    int* indices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double eigenValue = 0.;


    for (int i = 0; i < numMyElements_Row; ++i)
    {
        numEntries = 0;
        for (int j = 0; j < columnCount; ++j)
        {
            eigenValue = rEigenMatrix.coeff(i, j);
            if (eigenValue != 0)
            {
                ++numEntries;
                entriesVector.push_back(eigenValue);
                indicesVector.push_back(j);
            }
        }

        entries = &entriesVector[0];
        indices = &indicesVector[0];

//        convMatrix.InsertGlobalValues(myGlobalNodeIndices_Row[i], numEntries, entries, indices);
        convMatrix.InsertMyValues(i, numEntries, entries, indices);

        entriesVector.clear();
        indicesVector.clear();
        numEntries = 0;
    }

//    convMatrix.FillComplete();
    convMatrix.FillComplete(rColMap, rRowMap);



    return convMatrix;
}


Epetra_Vector ConversionTools::convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Epetra_Map rMap)
{
    int oldVectorSize = int(rEigenVector.rows());
    Epetra_Vector newVector(rMap);
//    int* myGlobalNodeIndices = rMap.MyGlobalElements();
    int numMyElements_Row = rMap.NumMyElements();


    double* vals = new double[1];
    int* indices = new int[1];
    for (int i = 0; i < oldVectorSize; ++i)
    {
        vals[0] = rEigenVector(i, 0);
//        indices[0] = myGlobalNodeIndices[i];
        indices[0] = i;
//        newVector.ReplaceGlobalValues(1, vals, indices);
        newVector.ReplaceMyValues(1, vals, indices);
    }

    delete[] vals;
    delete[] indices;

    return newVector;

}

