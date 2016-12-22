#include "base/Exception.h"
#include "metamodel/MetamodelException.h"
#include "math/FullMatrix.h"

#include "metamodel/PolynomialLeastSquaresFitting.h"

int main()
{
    try
    {

        int NSupportPoints = 21;

        NuTo::PolynomialLeastSquaresFitting PLSQF;

        NuTo::FullVector<double,Eigen::Dynamic> xVec;
        NuTo::FullVector<double,Eigen::Dynamic> yVec;

        NuTo::FullVector<double,Eigen::Dynamic> FixedCoeffs({2, 9,-2,7,-4,5});
        NuTo::FullVector<double,Eigen::Dynamic> CoeffGiven;
        NuTo::FullVector<double,Eigen::Dynamic> CoeffCalc;

        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> Result;


        xVec.resize(NSupportPoints);
        yVec.resize(NSupportPoints);

        for(int i=0; i<NSupportPoints; i++)
        {
            xVec(i) = double(i)/double((NSupportPoints-1));
            yVec(i) = double(i)/double((NSupportPoints-1));
        }


        // ---------------
        // test degree = 0
        // ---------------


        PLSQF.SetDegree(0);
        PLSQF.SetSupportPoints(1,1,xVec.Trans(),yVec.Trans());
        PLSQF.BuildDerived();

        CoeffCalc = PLSQF.GetPolynomialCoefficients();
        PLSQF.SolveTransformed(xVec.Trans(),Result);

        // Check coefficients
        if (    CoeffCalc.GetNumRows()>1 ||
                CoeffCalc.GetNumRows()<0 ||
                CoeffCalc(0) !=0.5)
        {
            throw NuTo::MetamodelException("[Testfile: PolynomialLeastSquaresFitting] Constant fit test: Calculated Coefficients are wrong");
        }
        //Check calues
        for(int i=0; i<Result.GetNumColumns();i++)
        {
            if (Result(i) != 0.5)
            {
                throw NuTo::MetamodelException("[Testfile: PolynomialLeastSquaresFitting] Constant fit test: Calculated results are wrong");
            }
        }


        // -----------------
        // test degree = 1-5
        // -----------------


        for(int pdegree=1; pdegree<6; pdegree++)
        {

            // Getting polynom coefficients
            CoeffGiven.resize(pdegree+1);
            for(int i=0; i<CoeffGiven.GetNumRows(); i++)
            {
                CoeffGiven(i) = FixedCoeffs(i);
            }

            // calculating function values
            for(int i=0; i<NSupportPoints;i++)
            {
                yVec(i) = 0;
                for(int j=0; j<CoeffGiven.GetNumRows(); j++)
                {
                    yVec(i) += CoeffGiven(j) * pow(xVec(i),j);
                }
            }

            PLSQF.SetDegree(pdegree);
            PLSQF.SetSupportPoints(1,1,xVec.Trans(),yVec.Trans());
            PLSQF.BuildDerived();

            CoeffCalc = PLSQF.GetPolynomialCoefficients();
            // Check coefficients
            if (    CoeffCalc.GetNumRows()>pdegree+1 ||
                    CoeffCalc.GetNumRows()<0)
            {
                throw NuTo::MetamodelException("[Testfile: PolynomialLeastSquaresFitting] test for degrees 1-5: Wrong number of coefficients");
            }
            // Check values
            for(int i=0; i<pdegree+1; i++)
            {
                if (std::abs(CoeffCalc(i)-CoeffGiven(i))>1e-6)
                {

                    throw NuTo::MetamodelException("[Testfile: PolynomialLeastSquaresFitting] test for degrees 1-5: Calculated polynomial coefficients differ to much from the given ones");
                }
            }


            PLSQF.SolveTransformed(xVec.Trans(),Result);


            for(int i=0; i<NSupportPoints; i++)
            {
                if(std::abs(Result(i) - yVec(i))>1e-6)
                {
                    throw NuTo::MetamodelException("[Testfile: PolynomialLeastSquaresFitting] test for degrees 1-5: Calculated function values differ to much from the given ones");
                }
            }
        }


        // ------------------------
        // Test boundary conditions
        // ------------------------


        PLSQF.AddBoundaryCondition(0,1.0);
        PLSQF.AddBoundaryCondition(0.5,1.0);
        PLSQF.AddBoundaryCondition(1.0,1.0);
        PLSQF.BuildDerived();
        xVec.resize(3);
        xVec(0) = 0;
        xVec(1) = 0.5;
        xVec(2) = 1.0;
        PLSQF.SolveTransformed(xVec.Trans(),Result);
        // Check values
        for(int i=0; i<3; i++)
        {
            if (std::abs(Result(i) - 1.0) > 1e-6)
            {
                throw NuTo::MetamodelException("[Testfile: PolynomialLeastSquaresFitting] test for boundary conditions: Calculated function values differ to much from set boundary conditions");
            }
        }
    }
    catch(NuTo::Exception e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return 1;
    }
}
