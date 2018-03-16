#include "PolynomialLeastSquaresFitting.h"
#include "base/Exception.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace NuTo;

void PolynomialLeastSquaresFitting::AddBoundaryCondition(std::pair<double, double> boundaryCondition)
{
    if (int(mBoundaryConditions.size()) < mDegree)
        mBoundaryConditions.push_back(boundaryCondition);
    else
        throw Exception(__PRETTY_FUNCTION__,
                        "The number of boundary conditions can't be higher than the degree of the polynomial.");
}

void PolynomialLeastSquaresFitting::AddBoundaryCondition(double x, double y)
{
    AddBoundaryCondition(std::pair<double, double>(x, y));
}

void PolynomialLeastSquaresFitting::BuildDerived()
{
    // check dimension of output
    int dimOutput = this->mSupportPoints.GetDimOutput();
    if (dimOutput != 1)
        throw Exception(__PRETTY_FUNCTION__, "Dimension of output must be 1.");

    if (mDegree < 0)
        throw Exception(__PRETTY_FUNCTION__, "Degree of polynomial not set.");

    int N_SupportPoints = this->mSupportPoints.GetNumSupportPoints();

    mPolynomialCoeffs.resize(mDegree + 1);
    Eigen::SparseMatrix<double> lhs(mDegree + 1, mDegree + 1);
    Eigen::VectorXd rhs(mDegree + 1);
    rhs.setZero();

    for (int k = 0; k < N_SupportPoints; k++)
    {
        for (int i = 0; i <= mDegree; i++)
        {
            for (int j = 0; j <= mDegree; j++)
            {
                if (i >= int(mBoundaryConditions.size()))
                    lhs.coeffRef(i, j) += pow(mSupportPoints.GetOrigSupportPointsInput()(0, k), i + j);
            }
            rhs[i] += pow(mSupportPoints.GetOrigSupportPointsInput()(0, k), i) *
                      mSupportPoints.GetOrigSupportPointsOutput()(0, k);
        }
    }

    for (int i = 0; i < int(mBoundaryConditions.size()); i++)
    {
        for (int j = 0; j <= mDegree; j++)
        {
            lhs.coeffRef(i, j) += pow(mBoundaryConditions[i].first, j);
        }
        rhs[i] = mBoundaryConditions[i].second;
    }

    lhs.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(lhs);
    solver.factorize(lhs);
    mPolynomialCoeffs = solver.solve(rhs);
}

double PolynomialLeastSquaresFitting::GetDegree() const
{
    return mDegree;
}

Eigen::VectorXd PolynomialLeastSquaresFitting::GetPolynomialCoefficients() const
{
    return mPolynomialCoeffs;
}

void PolynomialLeastSquaresFitting::SetDegree(int degree)
{
    if (degree < 0)
        throw Exception(__PRETTY_FUNCTION__, "Degree must be non negative.");
    mDegree = degree;
}

void PolynomialLeastSquaresFitting::SolveTransformed(const Eigen::MatrixXd& inputCoordinates,
                                                     Eigen::MatrixXd& outputCoordinates) const
{
    if (inputCoordinates.rows() != 1)
        throw Exception(__PRETTY_FUNCTION__, "Dimension of input (number of rows) has to be 1.");

    int numCoefficients = this->mPolynomialCoeffs.rows();
    if (numCoefficients != mDegree + 1)
        throw Exception(__PRETTY_FUNCTION__, "Invalid number of polynomial coefficients. Build model first.");

    // prepare output
    int numSamples = inputCoordinates.cols();
    outputCoordinates.resize(1, numSamples);
    outputCoordinates.setZero();

    // calculate output
    for (int k = 0; k < numSamples; k++)
    {
        for (int i = 0; i <= mDegree; i++)
        {
            // y2(k) = y2(k) + Coeff(i) * x2(k)^(i-1);
            outputCoordinates(0, k) += mPolynomialCoeffs[i] * pow(inputCoordinates(0, k), i);
        }
    }
}

void PolynomialLeastSquaresFitting::SetSupportPoints(int dimInput, int dimOutput, Eigen::MatrixXd inputCoordinates,
                                                     Eigen::MatrixXd outputCoordinates)
{
    if (dimInput != inputCoordinates.rows())
        throw Exception(__PRETTY_FUNCTION__, "Dimension of input must be equal to number of rows in the input matrix.");

    if (dimOutput != outputCoordinates.rows())
        throw Exception(__PRETTY_FUNCTION__,
                        "Dimension of output  must be equal to number of rows in the output matrix.");

    if (outputCoordinates.cols() != inputCoordinates.cols())
        throw Exception(__PRETTY_FUNCTION__,
                        "Number of samples (number of columns) in input and output matrix must be identical .");

    mSupportPoints.SetSupportPoints(inputCoordinates, outputCoordinates);
}
