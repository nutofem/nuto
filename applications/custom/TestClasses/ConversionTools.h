#include <map>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>


class ConversionTools
{
public:
#ifdef HAVE_MPI
    ConversionTools(Epetra_MpiComm rComm) : mComm(rComm){}
#else
    ConversionTools(Epetra_SerialComm rComm) : mComm(rComm){}
#endif
    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, bool rAsGlobal = false, bool rFillComplete = false);

    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_Map rRowMap, Epetra_Map rColMap, bool rFillComplete = false);

    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_CrsGraph rGraph, bool rFillComplete = false);


    Epetra_Vector convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, bool rAsGlobal = false);

    Epetra_Vector convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Epetra_Map rMap);

    std::vector<double> convertEpetraVector2StdVector(Epetra_Vector rVector);

    std::vector<double> convertEpetraMultiVector2StdVector(Epetra_MultiVector rMultiVector, int rVectorIndex = 0, bool rWrtMap = true);

private:
    int* map2Array_Int(std::map<int, int> rMap);

    std::map<int, int> invertMap_int(std::map<int, int> rMap);

    int defineNumLocalElements(int rNumGlobalElements);

    void generateDefaultGlobalMaps(Epetra_Map& rRowMap, Epetra_Map& rColumnMap, int rNumGlobalRows, int rNumGlobalColumns);

    void generateDefaultGlobalMap(Epetra_Map& rMap, int rNumGlobalElements);


protected:
#ifdef HAVE_MPI
    Epetra_MpiComm mComm;
#else
    Epetra_SerialComm mComm;
#endif

};//class
