#include "TrilinosUtils.h"

#include <mpi.h>

#include <Epetra_DataAccess.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>
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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

namespace TrilinosUtils
{
#ifdef HAVE_MPI
Epetra_Map createLinearMap(Epetra_MpiComm rComm, int rLocalNumberOfElements, int rIndexBase)
{
    return createLinearMap_KnownGlobal(rComm, -1, rLocalNumberOfElements, rIndexBase);
}

Epetra_Map createLinearMap_KnownGlobal(Epetra_MpiComm rComm, int rGlobalNumberOfElements, int rLocalNumberOfElements, int rIndexBase)
{
    Epetra_Map linearMap(rGlobalNumberOfElements, rLocalNumberOfElements, rIndexBase, rComm);
    return linearMap;
}


Epetra_Map createSpecificMap(Epetra_MpiComm rComm, int rLocalNumberOfElements, int rIndexBase, int* rSpecificMapping)
{
    return createSpecificMap_KnownGlobal(rComm, -1, rLocalNumberOfElements, rIndexBase, rSpecificMapping);
}

Epetra_Map createSpecificMap_KnownGlobal(Epetra_MpiComm rComm, int rGlobalNumberOfElements, int rLocalNumberOfElements, int rIndexBase, int *rSpecificMapping)
{
    Epetra_Map specificMap(-1, rLocalNumberOfElements, rSpecificMapping, rIndexBase, rComm);
    return specificMap;
}

#else
Epetra_Map createLinearMap(Epetra_SerialComm rComm, int rLocalNumberOfElements, int rIndexBase)
{
    Epetra_Map linearMap(-1, rLocalNumberOfElements, rIndexBase, rComm);
    return linearMap;
}

Epetra_Map createLinearMap_KnownGlobal(Epetra_SerialComm rComm, int rGlobalNumberOfElements, int rLocalNumberOfElements, int rIndexBase)
{
    Epetra_Map linearMap(rGlobalNumberOfElements, rLocalNumberOfElements, rIndexBase, rComm);
    return linearMap;
}


Epetra_Map createSpecificMap(Epetra_SerialComm rComm, int rLocalNumberOfElements, int rIndexBase, int* rSpecificMapping)
{
    Epetra_Map specificMap(-1, rLocalNumberOfElements, rSpecificMapping, rIndexBase, rComm);
    return specificMap;
}

Epetra_Map createSpecificMap_KnownGlobal(Epetra_SerialComm rComm, int rGlobalNumberOfElements, int rLocalNumberOfElements, int rIndexBase, int *rSpecificMapping)
{
    Epetra_Map specificMap(-1, rLocalNumberOfElements, rSpecificMapping, rIndexBase, rComm);
    return specificMap;
}
#endif

}
