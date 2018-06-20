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

////! @brief "Solver" for sparse linear algebra
struct CallBackSolver
{
    CallBackSolver(std::string solver, const Eigen::SparseMatrix<double>& C)
        : mSolver(solver)
        , mC(C)
    {
    }

    Eigen::VectorXd Solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b)
    {
        return mC * EigenSparseSolve(A, b, mSolver);
    }

    std::string mSolver;
    const Eigen::SparseMatrix<double>& mC;
};

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
    int DoStep(Eigen::VectorXd& start, QuasiStaticProblem& callBack, double newGlobalTime,
               std::string solverType = "EigenSparseLU");

    //! Writes the current time, the mean dof values and the sum of the residual into out, only for the given dof type
    //! and given dof numbers
    //! @param u ... the current value (e.g. displacements)
    //! @param out ... output stream
    //! @param dofType ... dof type
    //! @param dofNumbers ... dof numbers that are considered
    //! @param callBack ... call back class containing routins for function evaluation, derivative evaluation and
    //! consraints
    void WriteTimeDofResidual(Eigen::VectorXd& u, std::ostream& out, DofType dofType, std::vector<int> dofNumbers,
                              QuasiStaticProblem& callBack);

private:
};

} /* NuTo */
