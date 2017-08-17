//#define BOOST_TEST_MODULE Trilinos_Utils_Tests
//#include <boost/test/unit_test.hpp>
#include "BoostUnitTest.h"

//#include <mpi.h>

//#include <Epetra_DataAccess.h>
//#include <Epetra_CrsMatrix.h>
//#include <Epetra_Vector.h>
//#include <Epetra_Map.h>
//#include <Epetra_LinearProblem.h>
//#include <AztecOO.h>
//#include <Amesos.h>
//#include <Amesos_BaseSolver.h>
//#include <Amesos_ConfigDefs.h>
//#include <BelosLinearProblem.hpp>
//#include <BelosBlockGmresSolMgr.hpp>
//#include <BelosPseudoBlockGmresSolMgr.hpp>
//#include <BelosEpetraAdapter.hpp>
//#include <Ifpack.h>
//#include <Ifpack_AdditiveSchwarz.h>
////#include <Ifpack2_AdditiveSchwarz.hpp>
//#include <Teuchos_GlobalMPISession.hpp>
//#include <Teuchos_ParameterList.hpp>
//#include <Teuchos_RCP.hpp>


//#ifdef HAVE_MPI
//#include <Epetra_MpiComm.h>
//#else
//#include <Epetra_SerialComm.h>
//#endif



//BOOST_AUTO_TEST_CASE(EpetraMap_Creation)
//{
////    Teuchos::GlobalMPISession(&boost::unit_test::framework::master_test_suite().argc,
////                              &boost::unit_test::framework::master_test_suite().argv);
////#ifdef HAVE_MPI
////    Epetra_MpiComm comm(MPI_COMM_WORLD);
////#else
////    Epetra_SerialComm comm;
////#endif

////    int pcount = comm.NumProc();
////    int pid = comm.MyPID();
//    int pcount  = 4;
//    int pid = 1;
//    int numLocalElements = 5;
//    int numGlobalElements = pcount * numLocalElements;

//    int* specificMapping = new int[numLocalElements];
//    for (int i = 0; i < numLocalElements; ++i)
//    {
//        specificMapping[i] = pcount*i+pid;
//    }

////    Epetra_Map linearMap(-1, numLocalElements, 0, comm);
////    Epetra_Map specificMap(-1, numLocalElements, specificMapping, 0, comm);

//    int maxLID = numLocalElements-1;
//    int maxGID = numGlobalElements-1;
//    int myMaxGID = pcount*(numLocalElements-1)+pid;
//    BOOST_CHECK_EQUAL(maxGID, myMaxGID);
////    BOOST_CHECK_EQUAL(specificMap.MaxLID(), maxLID);
////    BOOST_CHECK_EQUAL(specificMap.MaxMyGID(), myMaxGID);
////    BOOST_CHECK_EQUAL(specificMap.MaxAllGID(), maxGID);


//}


BOOST_AUTO_TEST_CASE(AnotherOne)
{
    BOOST_CHECK_EQUAL(4,5);
}
