#pragma once

#include "mechanics/tools/TimeDependentProblem.h"
#include "mechanics/constraints/Constraints.h"
#include "boost/optional.hpp"

namespace NuTo
{
class QuasistaticSolver
{
public:
    //! Ctor
    //! @param equations system of equations including Gradient(), Hessian0() and UpdateHistory()
    //! @param dof dof type
    QuasistaticSolver(TimeDependentProblem& equations, DofType dof);

    //! @param constraints linear constraints
    //! @param numIndepententDofs number of independent dofs to build the constraint matrix
    void SetConstraints(Constraint::Constraints constraints, int numIndependentDofs);

    //! sets the global time required for evaluating the constraint right hand side
    //! @param globalTime global time
    void SetGlobalTime(double globalTime);

    //! builds the trial system where its residual contains forces equivialent to the applied constraints from time step
    //! t_n to t_n+1
    //! @param x independent dof values corresponding to globalTime
    //! @param globalTime t_n
    //! @param timeStep t_n+1 - t_n
    std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> TrialSystem(const Eigen::VectorXd& x, double globalTime,
                                                                        double timeStep);

    //! calculates and stores the history variables for the state x
    //! @param x independent dof values
    void UpdateHistory(const Eigen::VectorXd& x);

    //! evaluates the residual R(x), part of NuTo::NewtonRaphson::Problem
    //! @param x independent dof values
    Eigen::VectorXd Residual(const Eigen::VectorXd& x);

    //! evaluates the derivative dR/dx, part of NuTo::NewtonRaphson::Problem
    //! @param x independent dof values
    Eigen::SparseMatrix<double> Derivative(const Eigen::VectorXd& x);


    //! evaluates the norm of R, part of NuTo::NewtonRaphson::Problem
    //! @param residual residual vector
    double Norm(const Eigen::VectorXd& residual) const;


    //! prints values during the newton iterations, part of NuTo::NewtonRaphson::Problem
    //! @param residual residual vector
    void Info(int i, const Eigen::VectorXd& x, const Eigen::VectorXd& r) const;


    //! tolerance for Norm(R), public member because it is part of NuTo::NewtonRaphson::Problem
    double mTolerance = 1.e-10;

    //! shut up the Info() method
    void SetQuiet()
    {
        mQuiet = true;
    }

    //! Updates mProblem to time `newGlobalTime` and saves the new state mX upon convergence
    //! @param newGlobalTime new global time
    //! @param solverType solver type from NuTo::EigenSparseSolve(...)
    //! @return number of iterations required by the newton algorithm or boost::none upon failure to converge
    boost::optional<int> DoStep(double newGlobalTime, std::string solverType = "EigenSparseLU");

private:
    GlobalDofVector ToGlobalDofVector(const Eigen::VectorXd& x);

    Eigen::VectorXd mX;

    TimeDependentProblem& mProblem;
    Constraint::Constraints mConstraints;
    DofType mDof;
    Eigen::SparseMatrix<double> mCmat;

    double mGlobalTime = 0;
    double mTimeStep = 0;

    bool mQuiet = false;
};

} /* NuTo */
