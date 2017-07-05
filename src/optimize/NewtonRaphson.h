#pragma once

#include "optimize/NonlinearSolverBase.h"

namespace NuTo
{
//! @author Vitaliy Kindrachuk
//! @date December 2015
//! @brief ... standard class for solving system of nonlinear equations
class NewtonRaphson : public NonlinearSolverBase
{

public:
    //! @brief constructor
    NewtonRaphson();

    //! @brief sets the pointer to the analytic derivative of the residual function
    //! @param rParam ... parameters necessary to evaluate the residual
    //! @param rUnknown ... unknown vector
    void SetResidualDerivativeFunction(
            boost::function<Eigen::MatrixXd(const Eigen::VectorXd&, Eigen::VectorXd)> rResidualDerivativeFunction)
    {
        mResidualDerivativeFunctionBoost = rResidualDerivativeFunction;
        mAssignResidualResidualDerivative = true;
    }

    //! @brief returns mCheckNewtonRaphson, specifying iterator exit
    double GetCheckNewtonRaphson() const
    {
        return mCheckNewtonRaphson;
    }


    //! @brief perform iteration
    //! @param rUnknown ... unknown vector
    void Solve(Eigen::VectorXd& rUnknown) override;

    //! @brief ... the routine performs line search correction of the Newton step
    void LineSearch(Eigen::VectorXd& rXold, const double rFold, Eigen::VectorXd& rG, Eigen::VectorXd& rP,
                    Eigen::VectorXd& rX, double& rF, const double rStpmax, bool& rCheck, Eigen::VectorXd& rFvec) const;

    //! @brief ... the routine performs Newton-Raphson integration
    void NewtonRaphsonIterator(Eigen::VectorXd& x, bool& check) const;

protected:
    // pointer to the analytical derivative of the residual function (analytic Jacobi matrix)
    Eigen::MatrixXd (*mResidualDerivativeFunction)(const Eigen::VectorXd&, Eigen::VectorXd);

    boost::function<Eigen::MatrixXd(const Eigen::VectorXd&, Eigen::VectorXd)> mResidualDerivativeFunctionBoost;

    // a boolean variable which gets true if a jacobi function is assigned
    bool mAssignResidualResidualDerivative;

    //! @brief ... specify the exit: false on a normal exit; true if this is a local minimum of Fmin (that is the
    //! resudual function is not necessary zeroed)
    bool mCheckNewtonRaphson;
};
} // namespace NuTo
