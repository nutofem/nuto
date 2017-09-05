#include "base/Exception.h"
#include "metamodel/PolynomialLeastSquaresFitting.h"

int main()
{
    int NSupportPoints = 21;

    NuTo::PolynomialLeastSquaresFitting PLSQF;

    Eigen::VectorXd FixedCoeffs(6);
    FixedCoeffs << 2, 9, -2, 7, -4, 5;
    Eigen::VectorXd CoeffGiven;
    Eigen::VectorXd CoeffCalc;

    Eigen::MatrixXd Result;

    Eigen::VectorXd xVec(NSupportPoints);
    Eigen::VectorXd yVec(NSupportPoints);

    for (int i = 0; i < NSupportPoints; i++)
    {
        xVec(i) = double(i) / double((NSupportPoints - 1));
        yVec(i) = double(i) / double((NSupportPoints - 1));
    }


    // ---------------
    // test degree = 0
    // ---------------
    PLSQF.SetDegree(0);
    PLSQF.SetSupportPoints(1, 1, xVec.transpose(), yVec.transpose());
    PLSQF.BuildDerived();

    CoeffCalc = PLSQF.GetPolynomialCoefficients();
    PLSQF.SolveTransformed(xVec.transpose(), Result);

    // Check coefficients
    if (CoeffCalc.rows() > 1 || CoeffCalc.rows() < 0 || CoeffCalc(0) != 0.5)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Constant fit test: Calculated Coefficients are wrong");
    }
    // Check calues
    for (int i = 0; i < Result.cols(); i++)
    {
        if (Result(i) != 0.5)
        {
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Constant fit test: Calculated results are wrong");
        }
    }


    // -----------------
    // test degree = 1-5
    // -----------------
    for (int pdegree = 1; pdegree < 6; pdegree++)
    {
        // Getting polynom coefficients
        CoeffGiven.resize(pdegree + 1);
        for (int i = 0; i < CoeffGiven.rows(); i++)
        {
            CoeffGiven(i) = FixedCoeffs(i);
        }

        // calculating function values
        for (int i = 0; i < NSupportPoints; i++)
        {
            yVec(i) = 0;
            for (int j = 0; j < CoeffGiven.rows(); j++)
            {
                yVec(i) += CoeffGiven(j) * pow(xVec(i), j);
            }
        }

        PLSQF.SetDegree(pdegree);
        PLSQF.SetSupportPoints(1, 1, xVec.transpose(), yVec.transpose());
        PLSQF.BuildDerived();

        CoeffCalc = PLSQF.GetPolynomialCoefficients();
        // Check coefficients
        if (CoeffCalc.rows() > pdegree + 1 || CoeffCalc.rows() < 0)
        {
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Test for degrees 1-5: Wrong number of coefficients");
        }
        // Check values
        for (int i = 0; i < pdegree + 1; i++)
        {
            if (std::abs(CoeffCalc(i) - CoeffGiven(i)) > 1e-6)
            {
                throw NuTo::Exception(
                        __PRETTY_FUNCTION__,
                        "Test for degrees 1-5: Calculated polynomial coefficients differ to much from the given ones");
            }
        }

        PLSQF.SolveTransformed(xVec.transpose(), Result);
        for (int i = 0; i < NSupportPoints; i++)
        {
            if (std::abs(Result(i) - yVec(i)) > 1e-6)
            {
                throw NuTo::Exception(
                        __PRETTY_FUNCTION__,
                        "Test for degrees 1-5: Calculated function values differ to much from the given ones");
            }
        }
    }


    // ------------------------
    // Test boundary conditions
    // ------------------------
    PLSQF.AddBoundaryCondition(0, 1.0);
    PLSQF.AddBoundaryCondition(0.5, 1.0);
    PLSQF.AddBoundaryCondition(1.0, 1.0);
    PLSQF.BuildDerived();
    xVec.resize(3);
    xVec(0) = 0;
    xVec(1) = 0.5;
    xVec(2) = 1.0;
    PLSQF.SolveTransformed(xVec.transpose(), Result);
    // Check values
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(Result(i) - 1.0) > 1e-6)
        {
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Test for boundary conditions: Calculated function "
                                                       "values differ to much from set boundary conditions");
        }
    }
}
