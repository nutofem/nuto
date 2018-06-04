#pragma once

#include "nuto/mechanics/solver/Solve.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/constraints/ReducedSolutionSpace.h"
#include "nuto/math/EigenSparseSolve.h"
#include <iosfwd>

using namespace NuTo::Constraint;

namespace NuTo
{

//! @brief "Solver" for sparse linear algebra
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

class ImplicitCallBack
{
public:
    ImplicitCallBack(TimeDependentProblem& s, ReducedSolutionSpace& reducedSolutionSpaceOperator,
                     double tolerance = 1.e-10);

    //! evaluates the residual R(u), part of NuTo::NewtonRaphson::Problem
    //! @param u all dof values
    Eigen::VectorXd Residual(const Eigen::VectorXd& u);

    //! evaluates the derivative dR/dx, part of NuTo::NewtonRaphson::Problem
    //! @param u all dof values
    Eigen::SparseMatrix<double> Derivative(const Eigen::VectorXd& u);

    //! evaluates the residual R(u, t ,dt), part of NuTo::QuasistaticSolver::TrialState
    //! @param u all dof values
    Eigen::VectorXd Residual(const Eigen::VectorXd& u, double t, double dt);

    //! evaluates the residual D(u, t ,dt), part of NuTo::QuasistaticSolver::TrialState
    //! @param u all dof values
    Eigen::SparseMatrix<double> Derivative(const Eigen::VectorXd& u, double t, double dt);

    //! evaluates the (C^T)*(f_full + K_full * deltaBrhs), part of NuTo::QuasistaticSolver::TrialState
    //! @param u all dof values
    Eigen::VectorXd TrialStateRHS(const Eigen::VectorXd& u, const Eigen::SparseMatrix<double>& K_full,
                                  Eigen::VectorXd& deltaBrhsEigen, double t, double dt);

    //! evaluates the norm of R, part of NuTo::NewtonRaphson::Problem
    //! @param residual residual vector
    double Norm(const Eigen::VectorXd& residual) const;

    //! calculates and stores the history variables for the state u
    //! @param u all dof values
    void UpdateHistory(const Eigen::VectorXd& u);

    //! prints values during the newton iterations, part of NuTo::NewtonRaphson::Problem
    //! @param r residual residual vector
    void Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const;

    //! sets the global time required for evaluating the constraint right hand side
    //! @param globalTime global time
    void SetGlobalTime(double globalTime);

    //! sets the operator for the mapping of ind dofs to all dofs etc. (C * u_ind = u_all)
    //! @param reducedSolutionSpaceOperator ... object of the class ReducedSolutionSpace
    void SetReducedSolutionSpaceOperator(ReducedSolutionSpace& reducedSolutionSpaceOperator);

    //    void FillDofVector(DofVector<double>& destination, const Eigen::VectorXd& source) const;

public:
    //! tolerance for Norm(R), public member because it is part of NuTo::NewtonRaphson::Problem
    double mTolerance = 1.e-10;

    double mGlobalTime = 0;
    double mTimeStep = 0;

    //! the problem containing routines for gradient and its derivatives
    TimeDependentProblem& mProblem;

    //! the class describing the tranformation by the constraint matrix to the reduced solution space
    ReducedSolutionSpace& mReducedSolutionSpaceOperator;
};

class QuasistaticSolver
{
public:
    QuasistaticSolver()
    {
    }

    //! sets the global time required for evaluating the constraint right hand side
    //! @param globalTime global time
    void SetGlobalTime(double globalTime);

    //! computes the trial state of the system
    //! @param newGlobalTime new time, for which the trial state is to be computed
    //! @param solver that allows to extract the constraint displacements from previous steps
    Eigen::VectorXd TrialState(Eigen::VectorXd& start, double newGlobalTime, ImplicitCallBack& callBack,
                               std::string solverType);

    //! Updates mProblem to time `newGlobalTime` and saves the new state mX upon convergence
    //! @param newGlobalTime new global time
    //! @param solverType solver type from NuTo::EigenSparseSolve(...)
    //! @return number of iterations required by the newton algorithm, throws upon failure to converge
    int DoStep(Eigen::VectorXd& start, ImplicitCallBack& callBack, double newGlobalTime,
               std::string solverType = "EigenSparseLU", double tolerance = 1.e-10);

    //! Writes the current time, the mean dof values and the sum of the residual into out, only for the given dof type
    //! and given dof numbers
    //! @param out output stream
    //! @param dofType dof type
    //! @param dofNumbers dof numbers that are considered
    void WriteTimeDofResidual(Eigen::VectorXd& u, std::ostream& out, DofType dofType, std::vector<int> dofNumbers,
                              ImplicitCallBack& callBack);

private:
};

} /* NuTo */
