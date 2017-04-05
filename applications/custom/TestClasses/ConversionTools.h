#include <map>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
//#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
//#include <Epetra_FEVector.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>


class ConversionTools
{
public:
    ConversionTools(){}

    int* map2Array_Int(std::map<int, int> rMap);

    std::map<int, int> invertMap_int(std::map<int, int> rMap);

    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_Map rRowMap, Epetra_Map rColMap);

    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_CrsGraph rGraph);

    Epetra_Vector convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Epetra_Map rMap);


private:

};//class
