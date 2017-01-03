
#include "math/FullVector.h"
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
void NuTo::PolynomialLeastSquaresFitting::AddBoundaryCondition(std::pair<double,double> rBoundaryCondition)
{
    if (int(mBoundaryConditions.size()) < mDegree)
    {
        mBoundaryConditions.push_back(rBoundaryCondition);
    }
    else
    {
        throw MetamodelException("[NuTo::PolynomialLeastSquaresFitting::AddBoundaryCondition] The number of boundary conditions can't be higher than the degree of the polynom.");
    }
}

//! @brief ... Adds a pair of x and y coordinates that should be matched by the polynom
//! @param rX ... the x value of the boundary Condition
//! @param rY ... the y value of the boundary Condition
void NuTo::PolynomialLeastSquaresFitting::AddBoundaryCondition(double rX, double rY)
{
    AddBoundaryCondition(std::pair<double, double>(rX,rY));
}

//! @brief determine regression parameters
void NuTo::PolynomialLeastSquaresFitting::BuildDerived()
{
    // check dimension of output
    int dimOutput = this->mSupportPoints.GetDimOutput();
    if(dimOutput != 1)
    {
        throw MetamodelException("[NuTo::PolynomialLeastSquaresFitting::BuildDerived] dimension of output must be 1.");
    }
    if(mDegree < 0)
    {
        throw MetamodelException("[NuTo::PolynomialLeastSquaresFitting::BuildDerived] Degree of polynom not set.");
    }
    try
    {
        int N_SupportPoints = this -> mSupportPoints.GetNumSupportPoints();

        mPolynomialCoeffs.Resize(mDegree+1);
        NuTo::SparseMatrixCSRGeneral<double> lhs(mDegree+1, mDegree+1);
        NuTo::FullVector <double, Eigen::Dynamic> rhs(mDegree+1);

        for(int k=0; k < N_SupportPoints; k++)
        {
            for(int i=0; i <= mDegree; i++)
            {
                for(int j=0; j <= mDegree; j++)
                {
                    if ( i >= int(mBoundaryConditions.size()))
                    {
                        lhs.AddValue(i,j, pow(mSupportPoints.GetOrigSupportPointsInput()(0,k),i+j));
                    }
                }
                rhs(i) += pow(mSupportPoints.GetOrigSupportPointsInput()(0,k),i) * mSupportPoints.GetOrigSupportPointsOutput()(0,k);
            }
        }
        for (int i=0; i< int(mBoundaryConditions.size()); i++)
        {
            for (int j=0; j<=mDegree; j++)
            {
                lhs.AddValue(i,j,pow(mBoundaryConditions[i].first, j));
            }
            rhs(i) = mBoundaryConditions[i].second;
        }
        lhs.SetOneBasedIndexing();


        NuTo::SparseDirectSolverMUMPS Solver;
        Solver.SetShowTime(false);

        Solver.Solve(lhs,rhs,mPolynomialCoeffs);

    }
    catch(NuTo::Exception& e)
    {
        NuTo::MetamodelException myException(e.ErrorMessage());
        myException.AddMessage("[NuTo::PolynomialLeastSquaresFitting::BuildDerived] error calculating the polynomial coefficients.");
        throw myException;
    }

}

//! @brief ... Gets the degree of the polynom
//! @return degree of the polynom
double NuTo::PolynomialLeastSquaresFitting::GetDegree() const
{
    return mDegree;
}

//! @brief ... Gets the calculated polynomial coefficients
//! @return ... polynomial coefficients
NuTo::FullVector<double, Eigen::Dynamic> NuTo::PolynomialLeastSquaresFitting::GetPolynomialCoefficients() const
{
    return mPolynomialCoeffs;
}

std::string NuTo::PolynomialLeastSquaresFitting::GetTypeId()const
{
    return std::string("PolynomialLeastSquaresFitting");
}


//! @brief ... Sets the degree of the polynom
//! @param rDegree degree of the polynom
void NuTo::PolynomialLeastSquaresFitting::SetDegree(int rDegree)
{
    if(rDegree < 0)
    {
        throw MetamodelException("[NuTo::PolynomialLeastSquaresFitting::SetDegree] Degree must be non negative.");
    }
    mDegree = rDegree;
}

//! @brief ... calculate approximation (in transformed space)
//! @param rInputCoordinates ... matrix of input data points (transformed)
//! @param rOutputCoordinates ... vector of output data (transformed)
void NuTo::PolynomialLeastSquaresFitting::SolveTransformed(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates)const
{
    if(rInputCoordinates.GetNumRows() != 1)
    {
        throw MetamodelException("[NuTo::PolynomialLeastSquaresFitting::SolveTransformed] Dimension of input (number of rows) has to be 1.");
    }
    int numCoefficients = this->mPolynomialCoeffs.GetNumRows();
    if(numCoefficients != mDegree + 1)
    {
        throw MetamodelException("[NuTo::PolynomialLeastSquaresFitting::SolveTransformed] invalid number of polynomial coefficients. Build model first.");
    }
    try
    {
        // prepare output
        int numSamples = rInputCoordinates.GetNumColumns();
        rOutputCoordinates.Resize(1, numSamples);

        // calculate output
        for (int k = 0; k < numSamples; k++)
        {
            for(int i = 0; i <= mDegree; i++)
            {
                rOutputCoordinates(0,k) += mPolynomialCoeffs[i] * pow(rInputCoordinates(0,k),i);
                //y2(k) = y2(k) + Coeff(i) * x2(k)^(i-1);
            }
        }

    }
    catch(NuTo::Exception& e)
    {
        NuTo::MetamodelException myException(e.ErrorMessage());
        myException.AddMessage("[NuTo::PolynomialLeastSquaresFitting::SolveTransformed] error calculating response.");
        throw myException;
    }
}
