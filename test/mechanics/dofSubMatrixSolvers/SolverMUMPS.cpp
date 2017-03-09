#include "SolveSystem.h"
#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"


BOOST_AUTO_TEST_CASE(SolverMUMPS)
{
    NuTo::SolverMUMPS s;
    SolveAndCheckSystem(s);
}


