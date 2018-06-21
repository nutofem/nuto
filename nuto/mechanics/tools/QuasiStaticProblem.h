#pragma once

#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/constraints/ReducedSolutionSpace.h"
#include <iosfwd>

using namespace NuTo::Constraint;

namespace NuTo
{

class QuasiStaticProblem
{
public:
    QuasiStaticProblem(TimeDependentProblem& s, ReducedSolutionSpace& reducedSolutionSpaceOperator,
                       double tolerance = 1.e-10);

    //! evaluates the residual R(u), part of NuTo::NewtonRaphson::Problem
    //! @param u ... independent dof values
    Eigen::VectorXd Residual(const Eigen::VectorXd& u);

    //! evaluates the full residual R(u), e.g. to evaluate reaction forces
    //! @param u ... independent dof values
    DofVector<double> FullResidual(const Eigen::VectorXd& u);

    //! evaluates the derivative dR/dx, part of NuTo::NewtonRaphson::Problem
    //! @param u ... all dof values
    Eigen::SparseMatrix<double> Derivative(const Eigen::VectorXd& u);

    //! computes the trial state of the system
    //! @param start ... starting vector
    //! @param newGlobalTime ... new time, for which the trial state is to be computed
    //! @param solverType ... string defining which solver is used
    //! @return computes the trial state (starting vector for the Newton iteration)
    Eigen::VectorXd TrialState(DofVector<double>& start, double newGlobalTime, std::string solverType);

    //! evaluates the norm of R, part of NuTo::NewtonRaphson::Problem
    //! @param residual ... residual vector
    double Norm(const Eigen::VectorXd& residual) const;

    //! calculates and stores the history variables for the state u
    //! @param u .. all dof values
    void UpdateHistory(const Eigen::VectorXd& u);

    //! prints values during the newton iterations, part of NuTo::NewtonRaphson::Problem
    //! @param i ... iteration #
    //! @param x ... value e.g. displacements
    //! @param r ... residual residual vector
    void Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const;

    //! sets the global time required for evaluating the constraint right hand side
    //! @param globalTime ... global time
    void SetGlobalTime(double globalTime);

    //! gets the global time
    //! @return globalTime ... global time
    double GetGlobalTime() const
    {
        return mGlobalTime;
    }

    //! sets the operator for the mapping of ind dofs to all dofs etc. (C * u_ind = u_all)
    //! @param reducedSolutionSpaceOperator ... object of the class ReducedSolutionSpace
    void SetReducedSolutionSpaceOperator(ReducedSolutionSpace& reducedSolutionSpaceOperator);

    //    void FillDofVector(DofVector<double>& destination, const Eigen::VectorXd& source) const;

    //! tolerance for Norm(R), public member because it is part of NuTo::NewtonRaphson::Problem
    //! I don't think this is optimal (instead use tol as a parameter when calling NewtonRaphson)
    double mTolerance = 1.e-10;

private:
    //! evaluates the residual R(u, t ,dt), part of NuTo::QuasistaticSolver::TrialState
    //! @param u all dof values
    //! @param t ... time
    //! @param dt ... time step
    Eigen::VectorXd Residual(const Eigen::VectorXd& u, double t, double dt);

    //! evaluates the residual D(u, t ,dt), part of NuTo::QuasistaticSolver::TrialState
    //! @param u ... all dof values
    //! @param t ... time
    //! @param dt ... time step
    Eigen::SparseMatrix<double> Derivative(const Eigen::VectorXd& u, double t, double dt);


    double mGlobalTime = 0;
    double mTimeStep = 0;

    //! the problem containing routines for gradient and its derivatives
    TimeDependentProblem& mProblem;

    //! the class describing the tranformation by the constraint matrix to the reduced solution space
    ReducedSolutionSpace& mReducedSolutionSpaceOperator;

    //! this includes the complete solution vector of all dofs (not only the ones we solve for)
    DofVector<double> mSolution;
};

} /* NuTo */
