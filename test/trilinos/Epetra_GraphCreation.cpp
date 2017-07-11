//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
#include "BoostUnitTest.h"

#include <mpi.h>

#include <Epetra_DataAccess.h>
#include <Epetra_Map.h>

#include "Utilities.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


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
    if (procID == 0)
    {
        std::cout << "[Test_EpetraGraphCreation]: Linear graph with StaticProfile = false" << std::endl;
    }
    linearGraph_staticFalse.Print(std::cout);

    Epetra_CrsGraph linearGraph_staticFalse_augment(Epetra_DataAccess::Copy, linearMap, nonZerosPerRow, setStaticProfile);
    errCode = createGraph(nonZerosPerRow+1, linearGraph_staticFalse_augment);
    BOOST_CHECK_EQUAL(errCode, 0);
    linearGraph_staticFalse_augment.FillComplete();
    BOOST_CHECK_CLOSE(double(linearGraph_staticFalse_augment.NumMyNonzeros()), (nonZerosPerRow+1)*(nonZerosPerRow+2)/2, 0.001);
    if (procID == 0)
    {
        std::cout << "\n\n[Test_EpetraGraphCreation]: Augmented linear graph with StaticProfile = false" << std::endl;
    }
    linearGraph_staticFalse_augment.Print(std::cout);

    setStaticProfile = true;
    Epetra_CrsGraph linearGraph_staticTrue_augment(Epetra_DataAccess::Copy, linearMap, nonZerosPerRow, setStaticProfile);
    errCode = createGraph(nonZerosPerRow+1, linearGraph_staticTrue_augment);
    BOOST_CHECK_LT(errCode, 0); //if static profile is set to true, it has to occure some error
    linearGraph_staticTrue_augment.FillComplete();
    BOOST_CHECK_CLOSE(double(linearGraph_staticTrue_augment.NumMyNonzeros()), nonZerosPerRow*(nonZerosPerRow+1)/2, 0.001);
    if (procID == 0)
    {
        std::cout << "\n\n[Test_EpetraGraphCreation]: Augmented linear graph with StaticProfile = true" << std::endl;
    }
    linearGraph_staticTrue_augment.Print(std::cout);
}



BOOST_AUTO_TEST_CASE(EpetraGraphCreation)
{
    MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    checkEpetraGraphCreation(comm);

    MPI_Finalize();
}
































