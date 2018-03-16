#pragma once

#include <vector>
#include "SupportPoints.h"

namespace NuTo
{
class PolynomialLeastSquaresFitting
{
public:
    //! @brief Adds a pair of x and y coordinates that should be matched by the polynomial
    //! @param boundaryCondition a pair of x and y coordinates that should be matched by the polynomial
    void AddBoundaryCondition(std::pair<double, double> boundaryCondition);

    //! @brief Adds a pair of x and y coordinates that should be matched by the polynomial
    //! @param x the x value of the boundary condition
    //! @param y the y value of the boundary condition
    void AddBoundaryCondition(double x, double y);

    //! @brief determine regression parameters
    void BuildDerived();

    //! @brief Gets the degree of the polynomial
    //! @return degree of the polynomial
    double GetDegree() const;

    //! @brief Gets the calculated polynomial coefficients
    //! @return polynomial coefficients
    Eigen::VectorXd GetPolynomialCoefficients() const;

    //! @brief Sets the degree of the polynomial
    //! @param degree degree of the polynomial
    void SetDegree(int degree);

    //! @brief calculate approximation (in transformed space)
    //! @param inputCoordinates matrix of input data points (transformed)
    //! @param outputCoordinates vector of output data (transformed)
    void SolveTransformed(const Eigen::MatrixXd& inputCoordinates, Eigen::MatrixXd& outputCoordinates) const;

    void SetSupportPoints(int dimInput, int dimOutput, Eigen::MatrixXd inputCoordinates,
                          Eigen::MatrixXd outputCoordinates);

private:
    SupportPoints mSupportPoints;
    Eigen::VectorXd mPolynomialCoeffs;
    std::vector<std::pair<double, double>> mBoundaryConditions;
    int mDegree = -1;
};

} // namespace NuTo
