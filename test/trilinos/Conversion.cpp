//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
#include "BoostUnitTest.h"

#include <mpi.h>

#include <Epetra_DataAccess.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>

#include "Utilities.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#ifdef HAVE_MPI
void checkVectorConversions(Epetra_MpiComm rComm)
#else
void checkVectorConversions(Epetra_SerialComm rComm)
#endif
{
    int procID = rComm.MyPID();


    Eigen::Matrix<double, numLocalElements, 1> eigenVec;

    for (int i = 0; i < numLocalElements; ++i)
    {
        eigenVec(i,0) = i + procID*numLocalElements*10;
    }


    ConversionTools converter(rComm);

    Epetra_Map linearMap(-1, numLocalElements, 0, rComm);
    Epetra_Vector epetraVec = converter.convertEigen2EpetraVector(eigenVec, linearMap);

    epetraVec.Print(std::cout);

    BOOST_CHECK_EQUAL(epetraVec[numLocalElements-1], numLocalElements-1+procID*numLocalElements*10);
}


#ifdef HAVE_MPI
void checkMatrixConversions(Epetra_MpiComm rComm)
#else
void checkMatrixConversions(Epetra_SerialComm rComm)
#endif
{
    int errCode = 0;

    Epetra_Map linearMap(-1, numLocalElements, 0, rComm);
    Epetra_Map inverseLinearMap(-1, numLocalElements, 0, rComm);

    errCode = createInverseLinearMap(inverseLinearMap);
    BOOST_CHECK_EQUAL(errCode, 0);

    Eigen::SparseMatrix<double> eigenMatrix(numLocalElements, numLocalElements);
//    eigenMatrix.setIdentity();
    typedef Eigen::Triplet<double> doubleTrip;
    std::vector<doubleTrip> list;
    for (int j = 0; j < numLocalElements; ++j)
    {
        list.push_back(doubleTrip(j, j, 1));
    }
    list.push_back(doubleTrip(0, numLocalElements-1,3.141592));
    eigenMatrix.setFromTriplets(list.begin(), list.end());

    ConversionTools converter(rComm);
    Epetra_CrsMatrix epetraMatrix = converter.convertEigen2EpetraCrsMatrix(eigenMatrix, linearMap, inverseLinearMap);

    epetraMatrix.Print(std::cout);

    double* epetraValues = epetraMatrix[0];
    BOOST_CHECK_EQUAL(epetraValues[1], 3.141592);
}



BOOST_AUTO_TEST_CASE(Conversion)
{
    MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    checkVectorConversions(comm);

    checkMatrixConversions(comm);


    MPI_Finalize();
}
