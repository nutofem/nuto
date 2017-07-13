#include "metamodel/PolynomialLeastSquaresFitting.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::PolynomialLeastSquaresFitting::PolynomialLeastSquaresFitting()
    : Metamodel()
{
}

//! @brief ... Adds a pair of x and y coordinates that should be matched by the polynom
//! @param rBoundaryCondition ... a pair of x and y coordinates that should be matched by the polynom
void NuTo::PolynomialLeastSquaresFitting::AddBoundaryCondition(std::pair<double, double> rBoundaryCondition)
{
    if (int(mBoundaryConditions.size()) < mDegree)
    {
        mBoundaryConditions.push_back(rBoundaryCondition);
    }
    else
    {
        throw Exception("[NuTo::PolynomialLeastSquaresFitting::AddBoundaryCondition] The number of boundary "
                        "conditions can't be higher than the degree of the polynom.");
    }
}

//! @brief ... Adds a pair of x and y coordinates that should be matched by the polynom
//! @param rX ... the x value of the boundary Condition
//! @param rY ... the y value of the boundary Condition
void NuTo::PolynomialLeastSquaresFitting::AddBoundaryCondition(double rX, double rY)
{
    AddBoundaryCondition(std::pair<double, double>(rX, rY));
}

//! @brief determine regression parameters
void NuTo::PolynomialLeastSquaresFitting::BuildDerived()
{
    // check dimension of output
    int dimOutput = this->mSupportPoints.GetDimOutput();
    if (dimOutput != 1)
    {
        throw Exception("[NuTo::PolynomialLeastSquaresFitting::BuildDerived] dimension of output must be 1.");
    }
    if (mDegree < 0)
    {
        throw Exception("[NuTo::PolynomialLeastSquaresFitting::BuildDerived] Degree of polynom not set.");
    }

    int N_SupportPoints = this->mSupportPoints.GetNumSupportPoints();

    mPolynomialCoeffs.resize(mDegree + 1);
    NuTo::SparseMatrixCSRGeneral<double> lhs(mDegree + 1, mDegree + 1);
    Eigen::VectorXd rhs(mDegree + 1);
    rhs.setZero();

    for (int k = 0; k < N_SupportPoints; k++)
    {
        for (int i = 0; i <= mDegree; i++)
        {
            for (int j = 0; j <= mDegree; j++)
            {
                if (i >= int(mBoundaryConditions.size()))
                {
                    lhs.AddValue(i, j, pow(mSupportPoints.GetOrigSupportPointsInput()(0, k), i + j));
                }
            }
            rhs(i) += pow(mSupportPoints.GetOrigSupportPointsInput()(0, k), i) *
                      mSupportPoints.GetOrigSupportPointsOutput()(0, k);
        }
    }
    for (int i = 0; i < int(mBoundaryConditions.size()); i++)
    {
        for (int j = 0; j <= mDegree; j++)
        {
            lhs.AddValue(i, j, pow(mBoundaryConditions[i].first, j));
        }
        rhs(i) = mBoundaryConditions[i].second;
    }
    lhs.SetOneBasedIndexing();


    NuTo::SparseDirectSolverMUMPS Solver;
    Solver.SetShowTime(false);

    Solver.Solve(lhs, rhs, mPolynomialCoeffs);
}

//! @brief ... Gets the degree of the polynom
//! @return degree of the polynom
double NuTo::PolynomialLeastSquaresFitting::GetDegree() const
{
    return mDegree;
}

//! @brief ... Gets the calculated polynomial coefficients
//! @return ... polynomial coefficients
Eigen::VectorXd NuTo::PolynomialLeastSquaresFitting::GetPolynomialCoefficients() const
{
    return mPolynomialCoeffs;
}

//! @brief ... Sets the degree of the polynom
//! @param rDegree degree of the polynom
void NuTo::PolynomialLeastSquaresFitting::SetDegree(int rDegree)
{
    if (rDegree < 0)
    {
        throw Exception("[NuTo::PolynomialLeastSquaresFitting::SetDegree] Degree must be non negative.");
    }
    mDegree = rDegree;
}

//! @brief ... calculate approximation (in transformed space)
//! @param rInputCoordinates ... matrix of input data points (transformed)
//! @param rOutputCoordinates ... vector of output data (transformed)
void NuTo::PolynomialLeastSquaresFitting::SolveTransformed(const Eigen::MatrixXd& rInputCoordinates,
                                                           Eigen::MatrixXd& rOutputCoordinates) const
{
    if (rInputCoordinates.rows() != 1)
    {
        throw Exception("[NuTo::PolynomialLeastSquaresFitting::SolveTransformed] Dimension of input (number "
                        "of rows) has to be 1.");
    }
    int numCoefficients = this->mPolynomialCoeffs.rows();
    if (numCoefficients != mDegree + 1)
    {
        throw Exception("[NuTo::PolynomialLeastSquaresFitting::SolveTransformed] invalid number of polynomial "
                        "coefficients. Build model first.");
    }
    // prepare output
    int numSamples = rInputCoordinates.cols();
    rOutputCoordinates.resize(1, numSamples);
    rOutputCoordinates.setZero();

    // calculate output
    for (int k = 0; k < numSamples; k++)
    {
        for (int i = 0; i <= mDegree; i++)
        {
            rOutputCoordinates(0, k) += mPolynomialCoeffs[i] * pow(rInputCoordinates(0, k), i);
            // y2(k) = y2(k) + Coeff(i) * x2(k)^(i-1);
        }
    }
}
