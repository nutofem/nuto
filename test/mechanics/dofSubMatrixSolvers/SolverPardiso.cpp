#include "SolveSystem.h"
#include "mechanics/dofSubMatrixSolvers/SolverPardiso.h"


BOOST_AUTO_TEST_CASE(SolverPardiso)
{
    NuTo::SolverPardiso s(1);
    SolveAndCheckSystem(s);
}
