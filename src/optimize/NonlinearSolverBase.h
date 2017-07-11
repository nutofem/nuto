#pragma once

#include <boost/function.hpp>
#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @author Kindrachuk
//! @date December 2015
//! @brief ... standard abstract class for all solvers of nonlinear systems of equations
class NonlinearSolverBase
{
public:
    //! @brief constructor
    NonlinearSolverBase();

    //! @brief deconstructor
    virtual ~NonlinearSolverBase()
    {
    }

    //! @brief perform solving
    //! @param rUnknown ... unknown vector
    virtual void Solve(Eigen::VectorXd& rUnknown) = 0;

    //! @brief sets the pointer to the residual function
    //! @param rParam ... parameters necessary to evaluate the residual
    //! @param rUnknown ... unknown vector
    void
    SetResidualFunction(boost::function<Eigen::VectorXd(const Eigen::VectorXd&, Eigen::VectorXd)> rResidualFunction)
    {
        mResidualFunctionBoost = rResidualFunction;
        mAssignResidual = true;
    }

    //! @brief Sets the tolerance for the residual vector
    //! @remark The components of the residual should be of the same magnitude of order
    void SetResidualTolerance(double rTolResidual)
    {
        mTolResidual = rTolResidual;
    }

    //! @brief returns the tolerance for the residual vector
    double GetResidualTolerance() const
    {
        return mTolResidual;
    }

    //! @brief Sets the tolerance for the solution vector
    //! @remark The components of the rUnknown vector should be of the same magnitude of order
    void SetSolutionTolerance(double rTolSolution)
    {
        mTolSolution = rTolSolution;
    }

    //! @brief returns the tolerance for the solution rUnknown
    double GetSolutionTolerance() const
    {
        return mTolSolution;
    }

    //! @brief sets the  maximal number of iterations
    void SetMaxIterations(int rMaxIterations)
    {
        mMaxIterationsNumber = rMaxIterations;
    }

    //! @brief returns the tolerance for the residual vector
    double GetMaxIterations() const
    {
        return mMaxIterationsNumber;
    }

    //! @brief sets the list of parameters mParameter necessary to evaluate mResidualFunction
    void SetParameters(Eigen::VectorXd& rParameter)
    {
        mParameter = rParameter;
    }

    //! @brief ... numerical differentiation of the residual function mResidualFunction (numerical Jacobi matrix)
    //! @param rUnknown ... position at which the derivative is taken
    //! @param rFvec ... residual vector at the position rUnknown, rFvec = mResidualFunction(rParameter,rUnknown)
    Eigen::MatrixXd DResidualNum(Eigen::VectorXd rUnknown, Eigen::VectorXd& rFvec) const;

    //! @brief ... calculates 0.5*rFvec^2 and updates rFvec = mResidualFunction(mParameter,rUnknown)
    double Fmin(Eigen::VectorXd rUnknown, Eigen::VectorXd& rFvec) const;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info() const;

protected:
    // pointer to the residual function
    Eigen::VectorXd (*mResidualFunction)(const Eigen::VectorXd&, Eigen::VectorXd);

    boost::function<Eigen::VectorXd(const Eigen::VectorXd&, Eigen::VectorXd)> mResidualFunctionBoost;

    // tolerance for the residual vector
    double mTolResidual;

    // tolerance for the solution vector
    double mTolSolution;

    // maximal allowed number of iterations
    int mMaxIterationsNumber;

    // a boolean variable which gets true if a residual function is assigned
    bool mAssignResidual;

    // list of parameters necessary for evaluation of the resuduum
    Eigen::VectorXd mParameter;
};
} // namespace NuTo
