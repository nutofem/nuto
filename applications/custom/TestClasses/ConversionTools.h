#include <map>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
//#include <Tpetra_Distributor.hpp>

#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

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
    ConversionTools(Epetra_MpiComm rComm) : mComm(rComm)
    {
        mComm_teuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
//        mComm_teuchos = Teuchos::DefaultComm<int>::getComm();
//        mComm_teuchos = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }

    ConversionTools() : mComm(MPI_COMM_WORLD)
    {
        mComm_teuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
//        mComm_teuchos = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }
#else
    ConversionTools(Epetra_SerialComm rComm) : mComm(rComm){}
    ConversionTools(Teuchos::RCP<Teuchos::Comm<int> > rComm_teuchos) : mComm_teuchos(rComm_teuchos){}
#endif
#ifdef HAVE_MPI
    Epetra_MpiComm getComm_epetra(){return mComm;}
    Teuchos::RCP<const Teuchos::Comm<int>> getComm_tpetra(){return mComm_teuchos;}
#else
    Epetra_SerialComm getComm_epetra(){return mComm;}
#endif


    //******** Epetra conversion ********
    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, bool rAsGlobal = false, bool rFillComplete = false);

    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_Map rRowMap, Epetra_Map rColMap, bool rFillComplete = false);

    Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_CrsGraph rGraph, bool rFillComplete = false);

    Epetra_Vector convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, bool rAsGlobal = false);

    Epetra_Vector convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Epetra_Map rMap);

    std::vector<double> convertEpetraVector2StdVector(Epetra_Vector rVector);

    std::vector<double> convertEpetraMultiVector2StdVector(Epetra_MultiVector rMultiVector, int rVectorIndex = 0, bool rWrtMap = true);




    //******** Tpetra conversion ********
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> convertEigen2TpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, bool rAsGlobal = false, bool rFillComplete = false);

    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> convertEigen2TpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Teuchos::RCP<const Tpetra::Map<int, int>> rRowMap, Teuchos::RCP<const Tpetra::Map<int, int>> rColMap, bool rFillComplete = false);

    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> convertEigen2TpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Teuchos::RCP<const Tpetra::CrsGraph<int, int>> rGraph, bool rFillComplete = false);

    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int> > convertEigen2TpetraCrsMatrix_global(Eigen::SparseMatrix<double> rEigenMatrix, Teuchos::RCP<const Tpetra::Map<int, int> > rRowMap, Teuchos::RCP<const Tpetra::Map<int, int> > rColMap, bool rFillComplete);

    Teuchos::RCP<Tpetra::Vector<double, int, int>> convertEigen2TpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, bool rAsGlobal = false);

    Teuchos::RCP<Tpetra::Vector<double, int, int>> convertEigen2TpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Teuchos::RCP<const Tpetra::Map<int, int>> rMap);

    Teuchos::RCP<Tpetra::Vector<double, int, int> > convertEigen2TpetraVector_global(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Teuchos::RCP<const Tpetra::Map<int, int> > rMap);

    std::vector<double> convertTPetraVector2StdVector(Teuchos::RCP<const Tpetra::Vector<double, int, int>> rVector);

    std::vector<double> convertTpetraVector2StdVector(Tpetra::Vector<double, int, int> rVector);

    std::vector<double> convertTpetraMultiVector2StdVector(Teuchos::RCP<const Tpetra::MultiVector<double, int, int>> rMultiVector, int rVectorIndex = 0, bool rWrtMap = true);

    std::vector<double> convertTpetraMultiVector2StdVector(Tpetra::MultiVector<double, int, int> rMultiVector, int rVectorIndex = 0, bool rWrtMap = true);



    //******** Xpetra conversion ********
    template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
    Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>> convertEigen2XpetraCrsMatrix(Eigen::SparseMatrix<valuesType> rEigenMatrix, bool rAsGlobal = false, bool rFillComplete = false);

    template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
    Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>> convertEigen2XpetraCrsMatrix(Eigen::SparseMatrix<valuesType> rEigenMatrix, Teuchos::RCP<const Xpetra::Map<localOrdinalType, globalOrdinalType>> rRowMap, Teuchos::RCP<const Xpetra::Map<localOrdinalType, globalOrdinalType>> rColMap, bool rFillComplete = false);

    template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
    Teuchos::RCP<Xpetra::CrsMatrix<valuesType, localOrdinalType, globalOrdinalType>> convertEigen2XpetraCrsMatrix(Eigen::SparseMatrix<valuesType> rEigenMatrix, Teuchos::RCP<const Xpetra::CrsGraph<localOrdinalType, globalOrdinalType>> rGraph, bool rFillComplete = false);

    template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
    Teuchos::RCP<Xpetra::Vector<valuesType, localOrdinalType, globalOrdinalType>> convertEigen2XpetraVector(Eigen::Matrix<valuesType, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, bool rAsGlobal = false);

    template<typename valuesType, typename localOrdinalType, typename globalOrdinalType>
    Teuchos::RCP<Xpetra::Vector<valuesType, localOrdinalType, globalOrdinalType>> convertEigen2XpetraVector(Eigen::Matrix<valuesType, Eigen::Dynamic, Eigen::Dynamic> rEigenVector, Teuchos::RCP<const Xpetra::Map<localOrdinalType, globalOrdinalType>> rMap);

private:
    int* map2Array_Int(std::map<int, int> rMap);

    std::map<int, int> invertMap_int(std::map<int, int> rMap);

    int defineNumLocalElements(int rNumGlobalElements, int rRank, int rNumProcs);

    //******** Epetra support functions ********
    void generateDefaultGlobalMaps(Epetra_Map& rRowMap, Epetra_Map& rColumnMap, int rNumGlobalRows, int rNumGlobalColumns);

    void generateDefaultGlobalMap(Epetra_Map& rMap, int rNumGlobalElements);

    //******** Tpetra support functions ********
    void generateDefaultGlobalMapsTpetra(Teuchos::RCP<Tpetra::Map<int, int>>& rRowMap, Teuchos::RCP<Tpetra::Map<int, int>>& rColumnMap, int rNumGlobalRows, int rNumGlobalColumns);

    void generateDefaultGlobalMapsTpetra(Tpetra::Map<int, int>& rRowMap, Tpetra::Map<int, int>& rColumnMap, int rNumGlobalRows, int rNumGlobalColumns);

    void generateDefaultGlobalMapTpetra(Teuchos::RCP<Tpetra::Map<int, int>>& rMap, int rNumGlobalElements);

    void generateDefaultGlobalMapTpetra(Tpetra::Map<int, int>& rMap, int rNumGlobalElements);


    //******** Xpetra support functions ********
    template<typename localOrdinalType, typename globalOrdinalType>
    void generateDefaultGlobalMapsXpetra(Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>>& rRowMap, Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>>& rColumnMap, globalOrdinalType rNumGlobalRows, globalOrdinalType rNumGlobalColumns);

    template<typename localOrdinalType, typename globalOrdinalType>
    void generateDefaultGlobalMapXpetra(Teuchos::RCP<Xpetra::Map<localOrdinalType, globalOrdinalType>>& rMap, globalOrdinalType rNumGlobalElements);


protected:

    Teuchos::RCP<const Teuchos::Comm<int>> mComm_teuchos;
#ifdef HAVE_MPI
    Epetra_MpiComm mComm;
#else
    Epetra_SerialComm mComm;
#endif

};//class
