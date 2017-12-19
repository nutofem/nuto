#include "ConversionTools.h"
#include "PrintTools.h"

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
    int numMyElements_Row = rRowMap.NumMyElements();
    int* myGlobalIndices_Row = rRowMap.MyGlobalElements();
    int numMyElements_Col = rColMap.NumMyElements();
    int* myGlobalIndices_Col = rColMap.MyGlobalElements();



    Epetra_CrsMatrix convMatrix(Copy, rRowMap, rColMap, 0);


    std::vector<int> localIndicesVector;
    std::vector<int> globalIndicesVector;
    int* globalIndices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double value = 0.;

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

Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> ConversionTools::convertEigen2TpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, bool rAsGlobal, bool rFillComplete)
{
    int rowCount = rEigenMatrix.rows();
    int columnCount = rEigenMatrix.cols();

    if (rAsGlobal)
    {
        Teuchos::RCP<Tpetra::Map<int, int>> rowMap = Teuchos::rcp(new Tpetra::Map<int, int>(rowCount, 0, mComm_teuchos));
        Teuchos::RCP<Tpetra::Map<int, int>> columnMap = Teuchos::rcp(new Tpetra::Map<int, int>(columnCount, 0, mComm_teuchos));
        generateDefaultGlobalMapsTpetra(rowMap, columnMap, rowCount, columnCount);
        return convertEigen2TpetraCrsMatrix_global(rEigenMatrix, rowMap, columnMap, rFillComplete);
    }
    else
    {
        Teuchos::RCP<Tpetra::Map<int, int>> rowMap = Teuchos::rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), rowCount, 0, mComm_teuchos));
        Teuchos::RCP<Tpetra::Map<int, int>> columnMap = Teuchos::rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), columnCount, 0, mComm_teuchos));
        return convertEigen2TpetraCrsMatrix(rEigenMatrix, rowMap, columnMap, rFillComplete);
    }
}

Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> ConversionTools::convertEigen2TpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Teuchos::RCP<const Tpetra::Map<int, int> > rRowMap, Teuchos::RCP<const Tpetra::Map<int, int> > rColMap, bool rFillComplete)
{
    int numMyElements_Row = rEigenMatrix.rows();
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices_Row = rRowMap->getNodeElementList(); //rRowMap->getMyGlobalIndices();
    int numMyElements_Col = rEigenMatrix.cols();
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices_Col = rColMap->getNodeElementList();
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> convMatrix = Teuchos::rcp(new Tpetra::CrsMatrix<double, int, int>(rRowMap, rColMap, 0));


    std::vector<int> localIndicesVector;
    std::vector<int> globalIndicesVector;
    int* globalIndices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double value = 0.;


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

        convMatrix->insertGlobalValues(myGlobalIndices_Row[i], numEntries, entries, globalIndices);
//        convMatrix.InsertMyValues(i, numEntries, entries, indices);

        entriesVector.clear();
        localIndicesVector.clear();
        globalIndicesVector.clear();
        numEntries = 0;
    }


    if (rFillComplete)
        convMatrix->fillComplete(rColMap, rRowMap);

    return convMatrix;
}

Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> ConversionTools::convertEigen2TpetraCrsMatrix_global(Eigen::SparseMatrix<double> rEigenMatrix, Teuchos::RCP<const Tpetra::Map<int, int> > rRowMap, Teuchos::RCP<const Tpetra::Map<int, int> > rColMap, bool rFillComplete)
{
    int numMyElements_Row = rEigenMatrix.rows();
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices_Row = rRowMap->getNodeElementList(); //rRowMap->getMyGlobalIndices();
    int numMyElements_Col = rEigenMatrix.cols();
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices_Col = rColMap->getNodeElementList();
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> convMatrix = Teuchos::rcp(new Tpetra::CrsMatrix<double, int, int>(rRowMap, rColMap, 0));


    std::vector<int> localIndicesVector;
    std::vector<int> globalIndicesVector;
    int* globalIndices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double value = 0.;


    for (int i = 0; i < numMyElements_Row; ++i)
//    for (int i = 0; i < rowCount; ++i)
    {
        numEntries = 0;
        for (int j = 0; j < numMyElements_Col; ++j)
//        for (int j = 0; j < columnCount; ++j)
        {
            if (myGlobalIndices_Row[i] >= numMyElements_Row)
                std::cout << "WARNING: global ID [" << myGlobalIndices_Row[i] << "] too big for matrix row count " << numMyElements_Row << std::endl;
            else if (myGlobalIndices_Col[j] >= numMyElements_Col)
                std::cout << "WARNING: global ID [" << myGlobalIndices_Col[i] << "] too big for matrix column count " << numMyElements_Col << std::endl;
            else
            {
                value = rEigenMatrix.coeff(myGlobalIndices_Row[i], myGlobalIndices_Col[j]);
                if (value != 0)
                {
                    ++numEntries;
                    entriesVector.push_back(value);
    //                localIndicesVector.push_back(j);
                    globalIndicesVector.push_back(myGlobalIndices_Col[j]);
                }
            }
        }

        entries = &entriesVector[0];
//        localIndices = &localIndicesVector[0];
        globalIndices = &globalIndicesVector[0];

        convMatrix->insertGlobalValues(myGlobalIndices_Row[i], numEntries, entries, globalIndices);
//        convMatrix.InsertMyValues(i, numEntries, entries, indices);

        entriesVector.clear();
        localIndicesVector.clear();
        globalIndicesVector.clear();
        numEntries = 0;
    }


    if (rFillComplete)
        convMatrix->fillComplete(rColMap, rRowMap);

    return convMatrix;
}

Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> ConversionTools::convertEigen2TpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Teuchos::RCP<const Tpetra::CrsGraph<int, int> > rGraph, bool rFillComplete)
{
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> convertedMatrix = Teuchos::rcp(new Tpetra::CrsMatrix<double, int, int>(rGraph));
//    convertedMatrix->fillComplete();
    convertedMatrix->scale(0.0);

    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices_Row = rGraph->getMap()->getNodeElementList();
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices_Col = rGraph->getMap()->getNodeElementList();
    std::vector<int> globalIndicesVector;
    int* globalIndices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double value = 0.;

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

        convertedMatrix->sumIntoGlobalValues(myGlobalIndices_Row[i], numEntries, entries, globalIndices);

        entriesVector.clear();
        globalIndicesVector.clear();
        numEntries = 0;
    }

    if (rFillComplete)
        convertedMatrix->fillComplete();

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

Teuchos::RCP<Tpetra::Vector<double, int, int>> ConversionTools::convertEigen2TpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, bool rAsGlobal)
{
    int length = rEigenVector.rows();

    if (rAsGlobal)
    {
        Teuchos::RCP<Tpetra::Map<int, int>> map = Teuchos::rcp(new Tpetra::Map<int, int>(length, 0, mComm_teuchos));
        generateDefaultGlobalMapTpetra(map, length);

        return convertEigen2TpetraVector_global(rEigenVector, map);
    }
    else
    {
        Teuchos::RCP<const Tpetra::Map<int, int>> map = Teuchos::rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), length, 0, mComm_teuchos));

        return convertEigen2TpetraVector(rEigenVector, map);
    }
}

Teuchos::RCP<Tpetra::Vector<double, int, int>> ConversionTools::convertEigen2TpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Teuchos::RCP<const Tpetra::Map<int, int> > rMap)
{
    int vectorSize = int(rEigenVector.rows());
    Teuchos::RCP<Tpetra::Vector<double, int, int>> convertedVector = Teuchos::rcp(new Tpetra::Vector<double, int, int>(rMap));
    convertedVector->scale(0.0);
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices = rMap->getNodeElementList();

    for (int i = 0; i < vectorSize; ++i)
    {
//        convertedVector.SumIntoGlobalValues(1, values, globalIndices);
            convertedVector->replaceGlobalValue(myGlobalIndices[i], rEigenVector(i,0));
    }

//    delete myGlobalIndices;

    return convertedVector;
}

Teuchos::RCP<Tpetra::Vector<double, int, int>> ConversionTools::convertEigen2TpetraVector_global(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Teuchos::RCP<const Tpetra::Map<int, int> > rMap)
{
    int vectorSize = int(rEigenVector.rows());
    Teuchos::RCP<Tpetra::Vector<double, int, int>> convertedVector = Teuchos::rcp(new Tpetra::Vector<double, int, int>(rMap));
    convertedVector->scale(0.0);
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> myGlobalIndices = rMap->getNodeElementList();

    for (int i = 0; i < vectorSize; ++i)
    {
//        convertedVector.SumIntoGlobalValues(1, values, globalIndices);
        if (myGlobalIndices[i] >= vectorSize)
            std::cout << "WARNING: global ID [" << myGlobalIndices[i] << "] too big for vector of size " << vectorSize << std::endl;
        else
            convertedVector->replaceGlobalValue(myGlobalIndices[i], rEigenVector(myGlobalIndices[i],0));
    }

//    delete myGlobalIndices;

    return convertedVector;
}

std::vector<double> ConversionTools::convertTPetraVector2StdVector(Teuchos::RCP<const Tpetra::Vector<double, int, int>> rVector)
{
    return convertTpetraVector2StdVector(*rVector.get());
}

std::vector<double> ConversionTools::convertEpetraVector2StdVector(Epetra_Vector rVector)
{
    int n = rVector.MyLength();
    std::vector<double> convertedVector(n);

    for (int i = 0; i < n; ++i)
    {
        convertedVector[i] = rVector[i];
    }

    return convertedVector;
}

std::vector<double> ConversionTools::convertEpetraMultiVector2StdVector(Epetra_MultiVector rMultiVector, int rVectorIndex, bool rWrtMap)
{
    int local_n = rMultiVector.MyLength();
    int global_n = rMultiVector.GlobalLength();
    int* globalIDs = rMultiVector.Map().MyGlobalElements();
    std::vector<double> convertedVector;
    if (rWrtMap)
        convertedVector.resize(global_n);
    else
        convertedVector.resize(local_n);

    double* vals = rMultiVector[rVectorIndex];
    for (int i = 0; i < local_n; ++i)
    {
        if (rWrtMap)
            convertedVector[globalIDs[i]] = vals[i];
        else
            convertedVector[i] = vals[i];

    }

    return convertedVector;
}

std::vector<double> ConversionTools::convertTpetraVector2StdVector(Tpetra::Vector<double, int, int> rVector)
{
    int n = rVector.getLocalLength();
    std::vector<double> convertedVector(n);

    Teuchos::ArrayRCP<const double> vectorEntries = rVector.getData();

    for (int i = 0; i < n; ++i)
    {
        convertedVector[i] = vectorEntries[i];
    }

    return convertedVector;
}

std::vector<double> ConversionTools::convertTpetraMultiVector2StdVector(Teuchos::RCP<const Tpetra::MultiVector<double, int, int> > rMultiVector, int rVectorIndex, bool rWrtMap)
{
    return convertTpetraMultiVector2StdVector(*rMultiVector.get(), rVectorIndex, rWrtMap);
}

std::vector<double> ConversionTools::convertTpetraMultiVector2StdVector(Tpetra::MultiVector<double, int, int> rMultiVector, int rVectorIndex, bool rWrtMap)
{
    int local_n = rMultiVector.getLocalLength();
    int global_n = rMultiVector.getGlobalLength();
    Teuchos::ArrayView<const Tpetra::CrsMatrix<>::global_ordinal_type> globalIDs = rMultiVector.getMap()->getNodeElementList();
    std::vector<double> convertedVector;
    if (rWrtMap)
        convertedVector.resize(global_n);
    else
        convertedVector.resize(local_n);

    Teuchos::ArrayRCP<const double> vals = rMultiVector.getData(rVectorIndex);
    for (int i = 0; i < local_n; ++i)
    {
        if (rWrtMap)
            convertedVector[globalIDs[i]] = vals[i];
        else
            convertedVector[i] = vals[i];

    }

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

int ConversionTools::defineNumLocalElements(int rNumGlobalElements, int rRank, int rNumProcs)
{
    int numLocalElements = 0;

    double average = (double)rNumGlobalElements / rNumProcs;
    if (rRank < rNumProcs-1)
        numLocalElements = std::floor(average);
    else
        numLocalElements = rNumGlobalElements - (std::floor(average)*(rNumProcs-1));

    return numLocalElements;
}

void ConversionTools::generateDefaultGlobalMaps(Epetra_Map& rRowMap, Epetra_Map& rColumnMap, int rNumGlobalRows, int rNumGlobalColumns)
{
    int numLocalRows = defineNumLocalElements(rNumGlobalRows, mComm.MyPID(), mComm.NumProc());

    Epetra_Map newRowMap(rNumGlobalRows, numLocalRows, 0, mComm);
    Epetra_Map newColumnMap(rNumGlobalColumns, rNumGlobalColumns, 0, mComm);

    rRowMap = newRowMap;
    rColumnMap = newColumnMap;
}

void ConversionTools::generateDefaultGlobalMap(Epetra_Map& rMap, int rNumGlobalElements)
{
    int numLocalElements = defineNumLocalElements(rNumGlobalElements, mComm.MyPID(), mComm.NumProc());

    Epetra_Map newMap(rNumGlobalElements, numLocalElements, 0, mComm);

    rMap = newMap;
}

void ConversionTools::generateDefaultGlobalMapsTpetra(Teuchos::RCP<Tpetra::Map<int, int> > &rRowMap, Teuchos::RCP<Tpetra::Map<int, int> > &rColumnMap, int rNumGlobalRows, int rNumGlobalColumns)
{
    generateDefaultGlobalMapsTpetra(*rRowMap.get(), *rColumnMap.get(), rNumGlobalRows, rNumGlobalColumns);
}

void ConversionTools::generateDefaultGlobalMapsTpetra(Tpetra::Map<int, int> &rRowMap, Tpetra::Map<int, int> &rColumnMap, int rNumGlobalRows, int rNumGlobalColumns)
{
    int numLocalRows = defineNumLocalElements(rNumGlobalRows, mComm_teuchos->getRank(), mComm_teuchos->getSize());
    int numLocalCols = defineNumLocalElements(rNumGlobalColumns, mComm_teuchos->getRank(), mComm_teuchos->getSize());
    int* globalColIDs = new int[rNumGlobalColumns];
    for (int i = 0; i < rNumGlobalColumns; ++i)
        globalColIDs[i] = i;

    Tpetra::Map<int, int> newRowMap(rNumGlobalRows, numLocalRows, 0, mComm_teuchos);
    numLocalCols = rNumGlobalColumns;
    Tpetra::Map<int, int> newColumnMap(rNumGlobalColumns, numLocalCols, 0, mComm_teuchos);
//    Tpetra::Map<int, int> newColumnMap(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalColIDs, rNumGlobalColumns, 0, mComm_teuchos);

    rRowMap = newRowMap;
    rColumnMap = newColumnMap;
}

void ConversionTools::generateDefaultGlobalMapTpetra(Teuchos::RCP<Tpetra::Map<int, int> > &rMap, int rNumGlobalElements)
{
    generateDefaultGlobalMapTpetra(*rMap.get(), rNumGlobalElements);
}

void ConversionTools::generateDefaultGlobalMapTpetra(Tpetra::Map<int, int> &rMap, int rNumGlobalElements)
{
    int numLocalElements = defineNumLocalElements(rNumGlobalElements, mComm_teuchos->getRank(), mComm_teuchos->getSize());

    Tpetra::Map<int, int> newMap(rNumGlobalElements, numLocalElements, 0, mComm_teuchos);

    rMap = newMap;
}





template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType> > ConversionTools::convertEigen2XpetraCrsMatrix(Eigen::SparseMatrix<valuesType> rEigenMatrix, bool rAsGlobal, bool rFillComplete)
{
    int rowCount = rEigenMatrix.rows();
    int columnCount = rEigenMatrix.cols();

    if (rAsGlobal)
    {
        Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>> rowMap = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(rowCount, 0, mComm_teuchos));
        Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>> columnMap = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(columnCount, 0, mComm_teuchos));
        generateDefaultGlobalMapsXpetra(rowMap, columnMap, rowCount, columnCount);
        return convertEigen2XpetraCrsMatrix(rEigenMatrix, rowMap, columnMap, rFillComplete);
    }
    else
    {
        Teuchos::RCP<const Tpetra::Map<localOrdinalType, globalOrdinalType>> rowMap = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), rowCount, 0, mComm_teuchos));
        Teuchos::RCP<const Tpetra::Map<localOrdinalType, globalOrdinalType>> columnMap = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), columnCount, 0, mComm_teuchos));
        return convertEigen2XpetraCrsMatrix(rEigenMatrix, rowMap, columnMap, rFillComplete);
    }
}


template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>> ConversionTools::convertEigen2XpetraCrsMatrix(Eigen::SparseMatrix<valuesType> rEigenMatrix, Teuchos::RCP<const Xpetra::Map<localOrdinalType, globalOrdinalType>> rRowMap, Teuchos::RCP<const Xpetra::Map<localOrdinalType, globalOrdinalType>> rColMap, bool rFillComplete)
{
    int numMyElements_Row = rRowMap->getNodeNumElements();
    Teuchos::ArrayView<const globalOrdinalType> myGlobalIndices_Row = rRowMap->getNodeElementList(); //rRowMap->getMyGlobalIndices();
    int numMyElements_Col = rColMap->getNodeNumElements();
    Teuchos::ArrayView<const globalOrdinalType> myGlobalIndices_Col = rColMap->getNodeElementList();
    Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>> convMatrix = Teuchos::rcp(new Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>(rRowMap, rColMap, 0));


    std::vector<localOrdinalType> localIndicesVector;
    std::vector<globalOrdinalType> globalIndicesVector;
    int numEntries = 0;
    std::vector<valuesType> entriesVector;
    double value = 0.;

    for (int i = 0; i < numMyElements_Row; ++i)
    {
        numEntries = 0;
        for (int j = 0; j < numMyElements_Col; ++j)
        {
            value = rEigenMatrix.coeff(i, j);
            if (value != 0)
            {
                ++numEntries;
                entriesVector.push_back(value);
                globalIndicesVector.push_back(myGlobalIndices_Col[j]);
            }
        }

        Teuchos::ArrayView<globalOrdinalType> colInds = Teuchos::arrayViewFromVector(globalIndicesVector);
        Teuchos::ArrayView<valuesType> values = Teuchos::arrayViewFromVector(entriesVector);

        convMatrix->insertGlobalValues(myGlobalIndices_Row[i], colInds, values);

        entriesVector.clear();
        localIndicesVector.clear();
        globalIndicesVector.clear();
        numEntries = 0;
    }


    if (rFillComplete)
        convMatrix->fillComplete();

    return convMatrix;
}




template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>> ConversionTools::convertEigen2XpetraCrsMatrix(Eigen::SparseMatrix<valuesType> rEigenMatrix, Teuchos::RCP<const Xpetra::CrsGraph<localOrdinalType, globalOrdinalType>> rGraph, bool rFillComplete)
{
    Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>> convertedMatrix = Teuchos::rcp(new Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>(rGraph));
    convertedMatrix->fillComplete();
    convertedMatrix->scale(0.0);

    Teuchos::ArrayView<globalOrdinalType> myGlobalIndices_Row = rGraph->getMap()->getNodeElementList();
    Teuchos::ArrayView<globalOrdinalType> myGlobalIndices_Col = rGraph->getMap()->getNodeElementList();
    std::vector<globalOrdinalType> globalIndicesVector;
    int numEntries = 0;
    std::vector<valuesType> entriesVector;
    double value = 0.;

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


        Teuchos::ArrayView<globalOrdinalType> colInds = Teuchos::arrayViewFromVector(globalIndicesVector);
        Teuchos::ArrayView<valuesType> values = Teuchos::arrayViewFromVector(entriesVector);

//        convertedMatrix->sumIntoGlobalValues(myGlobalIndices_Row[i], numEntries, entries, globalIndices);
        convertedMatrix->insertGlobalValues(myGlobalIndices_Row[i], colInds, values);

        entriesVector.clear();
        globalIndicesVector.clear();
        numEntries = 0;
    }

    if (rFillComplete)
        convertedMatrix->fillComplete();

    return convertedMatrix;
}


template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
Teuchos::RCP<Xpetra::Vector<valuesType, localOrdinalType, globalOrdinalType>> ConversionTools::convertEigen2XpetraVector(Eigen::Matrix<valuesType, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Teuchos::RCP<const Xpetra::Map<localOrdinalType, globalOrdinalType>> rMap)
{
    int vectorSize = int(rEigenVector.rows());
    Teuchos::RCP<Xpetra::Vector<valuesType, localOrdinalType, globalOrdinalType>> convertedVector = Teuchos::rcp(new Xpetra::Vector<valuesType, localOrdinalType, globalOrdinalType>(rMap));
    convertedVector->scale(0.0);
    Teuchos::ArrayView<globalOrdinalType> myGlobalIndices = rMap->getNodeElementList();

    for (int i = 0; i < vectorSize; ++i)
    {
        convertedVector->replaceGlobalValue(myGlobalIndices[i], rEigenVector(i,0));
    }
    return convertedVector;
}

template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
Teuchos::RCP<Xpetra::Vector<valuesType, localOrdinalType, globalOrdinalType> > ConversionTools::convertEigen2XpetraVector(Eigen::Matrix<valuesType, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, bool rAsGlobal)
{
    int length = rEigenVector.rows();

    if (rAsGlobal)
    {
        Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>> map = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(length, 0, mComm_teuchos));
        generateDefaultGlobalMapTpetra(map, length);

        return convertEigen2XpetraVector(rEigenVector, map);
    }
    else
    {
        Teuchos::RCP<const Xpetra::Map<localOrdinalType, globalOrdinalType>> map = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), length, 0, mComm_teuchos));

        return convertEigen2XpetraVector(rEigenVector, map);
    }
}

template<typename localOrdinalType, typename globalOrdinalType>
void ConversionTools::generateDefaultGlobalMapsXpetra(Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType> > &rRowMap, Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType> > &rColumnMap, globalOrdinalType rNumGlobalRows, globalOrdinalType rNumGlobalColumns)
{
    int numLocalRows = defineNumLocalElements(rNumGlobalRows, mComm_teuchos->getRank(), mComm_teuchos->getSize());

    Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>> newRowMap = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(rNumGlobalRows, numLocalRows, 0, mComm_teuchos));
    Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>> newColumnMap = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(rNumGlobalColumns, rNumGlobalColumns, 0, mComm_teuchos));

    rRowMap = newRowMap;
    rColumnMap = newColumnMap;
}

template<typename localOrdinalType, typename globalOrdinalType>
void ConversionTools::generateDefaultGlobalMapXpetra(Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType> > &rMap, globalOrdinalType rNumGlobalElements)
{
    int numLocalElements = defineNumLocalElements(rNumGlobalElements, mComm_teuchos->getRank(), mComm_teuchos->getSize());

    Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>> newMap = Teuchos::rcp(new Xpetra::Map<localOrdinalType, globalOrdinalType>(rNumGlobalElements, numLocalElements, 0, mComm_teuchos));

    rMap = newMap;
}











