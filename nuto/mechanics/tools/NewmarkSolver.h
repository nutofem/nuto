#pragma once

#include "nuto/mechanics/solver/Solve.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include <iosfwd>

namespace NuTo
{
class NewmarkSolver
{
public:
    //! Ctor
    //! @param equations system of equations including Gradient(), Hessian0() and UpdateHistory()
    //! @param dof dof type
    NewmarkSolver(TimeDependentProblem<2>& equations, DofType dof, double beta = 0.25, double gamma = 0.5);

    //! Ctor
    //! @param equations system of equations including Gradient(), Hessian0() and UpdateHistory()
    //! @param dofs multiple dof types
    NewmarkSolver(TimeDependentProblem<2>& equations, std::vector<DofType> dofs, double beta = 0.25,
                  double gamma = 0.5);

    //! @param constraints linear constraints
    void SetConstraints(Constraint::Constraints constraints);

    //! sets the global time required for evaluating the constraint right hand side
    //! @param globalTime global time
    void SetGlobalTime(double globalTime);

    Eigen::SparseMatrix<double> Hessian(const Eigen::SparseMatrix<double>& hessian0,
                                        const Eigen::SparseMatrix<double>& hessian1,
                                        const Eigen::SparseMatrix<double>& hessian2, double delta_t) const;

    //! computes the trial state of the system
    //! @param newGlobalTime new time, for which the trial state is to be computed
    //! @param solver that allows to extract the constraint displacements from previous steps
    DofVector<double> TrialState(double newGlobalTime, const NuTo::ConstrainedSystemSolver& solver);

    //! calculates and stores the history variables for the state x
    //! @param x independent dof values
    void UpdateHistory(const DofVector<double>& x);

    //! evaluates the residual R(u), part of NuTo::NewtonRaphson::Problem
    //! @param u independent dof values
    DofVector<double> Residual(const DofVector<double>& u);

    //! evaluates the derivative dR/dx, part of NuTo::NewtonRaphson::Problem
    //! @param u independent dof values
    DofMatrixSparse<double> Derivative(const DofVector<double>& u);


    //! evaluates the norm of R, part of NuTo::NewtonRaphson::Problem
    //! @param residual residual vector
    double Norm(const DofVector<double>& residual) const;


    //! prints values during the newton iterations, part of NuTo::NewtonRaphson::Problem
    //! @param r residual residual vector
    void Info(int i, const DofVector<double>& x, const DofVector<double>& r) const;


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
    //! @return number of iterations required by the newton algorithm, throws upon failure to converge
    int DoStep(double newGlobalTime, std::string solverType = "EigenSparseLU");

    //! Writes the current time, the mean dof values and the sum of the residual into out, only for the given dof type
    //! and given dof numbers
    //! @param out output stream
    //! @param dofType dof type
    //! @param dofNumbers dof numbers that are considered
    void WriteTimeDofResidual(std::ostream& out, DofType dofType, std::vector<int> dofNumbers);

private:
    //! @var mX last updated dof state
    std::array<DofVector<double>, 3> mX;

    TimeDependentProblem<2>& mProblem;
    Constraint::Constraints mConstraints;

    std::vector<DofType> mDofs;
    DofMatrixSparse<double> mCmatUnit;

    double mGlobalTime = 0;
    double mTimeStep = 0;

    double mBeta;
    double mGamma;

    bool mQuiet = false;
};

} /* NuTo */
