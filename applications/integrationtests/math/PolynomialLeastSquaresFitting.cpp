#include "base/Exception.h"
#include "math/PolynomialLeastSquaresFitting.h"

using namespace NuTo;

int main()
{
    int numSupportPoints = 21;
    Eigen::VectorXd xVec(numSupportPoints);
    Eigen::VectorXd yVec(numSupportPoints);

    for (int i = 0; i < numSupportPoints; i++)
    {
        xVec(i) = double(i) / double((numSupportPoints - 1));
        yVec(i) = double(i) / double((numSupportPoints - 1));
    }

    // ---------------
    // test degree = 0
    // ---------------
    PolynomialLeastSquaresFitting PLSQF;
    PLSQF.SetDegree(0);
    PLSQF.SetSupportPoints(1, 1, xVec.transpose(), yVec.transpose());
    PLSQF.BuildDerived();

    auto coeffCalc = PLSQF.GetPolynomialCoefficients();
    Eigen::MatrixXd result;
    PLSQF.SolveTransformed(xVec.transpose(), result);

    // Check coefficients
    if (coeffCalc.rows() > 1 || coeffCalc.rows() < 0 || std::abs(coeffCalc[0] - 0.5) > 1.e-10)
        throw Exception(__PRETTY_FUNCTION__, "Constant fit test: Calculated Coefficients are wrong");

    // Check calues
    for (int i = 0; i < result.cols(); i++)
    {
        if (std::abs(result(i) - 0.5) > 1.e-10)
            throw Exception(__PRETTY_FUNCTION__, "Constant fit test: Calculated results are wrong");
    }

    // -----------------
    // test degree = 1-5
    // -----------------
    Eigen::VectorXd coeffGiven;
    Eigen::VectorXd fixedCoeffs(6);
    fixedCoeffs << 2, 9, -2, 7, -4, 5;
    for (int pdegree = 1; pdegree < 6; pdegree++)
    {
        // Getting polynom coefficients
        coeffGiven.resize(pdegree + 1);
        for (int i = 0; i < coeffGiven.rows(); i++)
        {
            coeffGiven(i) = fixedCoeffs(i);
        }

        // calculating function values
        for (int i = 0; i < numSupportPoints; i++)
        {
            yVec(i) = 0;
            for (int j = 0; j < coeffGiven.rows(); j++)
            {
                yVec(i) += coeffGiven(j) * pow(xVec(i), j);
            }
        }

        PLSQF.SetDegree(pdegree);
        PLSQF.SetSupportPoints(1, 1, xVec.transpose(), yVec.transpose());
        PLSQF.BuildDerived();

        coeffCalc = PLSQF.GetPolynomialCoefficients();
        // Check coefficients
        if (coeffCalc.rows() > pdegree + 1 || coeffCalc.rows() < 0)
            throw Exception(__PRETTY_FUNCTION__, "Test for degrees 1-5: Wrong number of coefficients");

        // Check values
        for (int i = 0; i < pdegree + 1; i++)
        {
            if (std::abs(coeffCalc(i) - coeffGiven(i)) > 1e-6)
                throw Exception(
                        __PRETTY_FUNCTION__,
                        "Test for degrees 1-5: Calculated polynomial coefficients differ too much from the given ones");
        }

        PLSQF.SolveTransformed(xVec.transpose(), result);
        for (int i = 0; i < numSupportPoints; i++)
        {
            if (std::abs(result(i) - yVec(i)) > 1e-6)
            {
                throw Exception(__PRETTY_FUNCTION__,
                                "Test for degrees 1-5: Calculated function values differ too much from the given ones");
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
    PLSQF.SolveTransformed(xVec.transpose(), result);

    // Check values
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(result(i) - 1.0) > 1e-6)
            throw Exception(__PRETTY_FUNCTION__, "Test for boundary conditions: Calculated function "
                                                 "values differ to much from set boundary conditions");
    }
}
