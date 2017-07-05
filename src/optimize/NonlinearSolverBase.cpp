#include "optimize/NonlinearSolverBase.h"
#include "optimize/OptimizeException.h"

using namespace std;

NuTo::NonlinearSolverBase::NonlinearSolverBase()
{
    mTolResidual = 1.0e-12;
    mTolSolution = numeric_limits<double>::epsilon();
    mMaxIterationsNumber = 100;
    mResidualFunction = nullptr;
    mAssignResidual = false;
}

Eigen::MatrixXd NuTo::NonlinearSolverBase::DResidualNum(Eigen::VectorXd rUnknown, Eigen::VectorXd& rFvec) const
{
    const double EPS = 1.0e-8;
    int n = rUnknown.rows();
    Eigen::MatrixXd deriv(n, n);
    Eigen::VectorXd xh = rUnknown;

    if (mResidualFunction == nullptr && mAssignResidual == false)
    {
        throw OptimizeException(__PRETTY_FUNCTION__, "The pointer to the residual function is required.");
    }

    for (int j = 0; j < n; j++)
    {
        double temp = xh[j];
        double h = EPS * std::abs(temp);
        if (h == 0.0)
            h = EPS;
        xh[j] = temp + h;
        h = xh[j] - temp;
        Eigen::VectorXd f = (mResidualFunctionBoost)(this->mParameter, xh);
        xh[j] = temp;
        for (int i = 0; i < n; i++)
            deriv(i, j) = (f[i] - rFvec[i]) / h;
    }
    return deriv;
}

double NuTo::NonlinearSolverBase::Fmin(Eigen::VectorXd rUnknown, Eigen::VectorXd& rFvec) const
{
    if (mResidualFunction == nullptr && mAssignResidual == false)
    {
        throw OptimizeException(__PRETTY_FUNCTION__, "The pointer to the residual function is required.");
    }

    rFvec = (mResidualFunctionBoost)(this->mParameter, rUnknown);
    double sum = 0;
    sum = rFvec.dot(rFvec);
    return 0.5 * sum;
}

void NuTo::NonlinearSolverBase::Info() const {}

