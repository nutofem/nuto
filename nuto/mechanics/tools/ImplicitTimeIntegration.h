#pragma once

#include "nuto/mechanics/solver/Solve.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/tools/QuasiStaticProblem.h"
#include "nuto/mechanics/constraints/ReducedSolutionSpace.h"
#include "nuto/math/EigenSparseSolve.h"
#include <iosfwd>

using namespace NuTo::Constraint;

namespace NuTo
{

class ImplicitTimeIntegration
{
public:
    ImplicitTimeIntegration()
    {
    }

    //! sets the global time required for evaluating the constraint right hand side
    //! @param globalTime global time
    void SetGlobalTime(double globalTime);

    //! Updates mProblem to time `newGlobalTime` and saves the new state mX upon convergence
    //! @param start ... starting vector
    //! @param implicitProblem ... call back class containing routins for function evaluation, derivative evaluation and
    //! constraints
    //! @param newGlobalTime ... new time, for which the trial state is to be computed
    //! @param solverType solver type from NuTo::EigenSparseSolve(...)
    //! @return number of iterations required by the newton algorithm, throws upon failure to converge
    int DoStep(DofVector<double>& start, QuasiStaticProblem& callBack, double newGlobalTime,
               std::string solverType = "EigenSparseLU");
};

} /* NuTo */
