// $Id$
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include "math/fortran_routines.h"
#include "math/MathException.h"
#include "math/NuToMath.h"

void NuTo::MatrixOperations::mult(const double *rA, int rNumRowsA, int rNumColumnsA,
                                  const double *rB, int rNumRowsB, int rNumColumnsB,
                                  char rTransA, char rTransB, double rAlpha,
                                  double *rResultC, double rBeta )
{
#ifdef HAVE_LAPACK
// just a simple wrapper file for the fortran dgemm to make the source code better readable
// for performance relevant parts, replace the mult routine directly with the dgemm routine
    try
    {
        if (rTransA=='n')
        {
            if (rTransB=='n')
            {
                //A*B
                if (rNumColumnsA!=rNumRowsB)
                    throw NuTo::MathException("[FullMatrix::mult] number of columns of A must be equal to number of rows of B.");
                dgemm_ (&rTransA, &rTransB, &rNumRowsA, &rNumColumnsB,  &rNumColumnsA, &rAlpha,
                        rA, &rNumRowsA, rB, &rNumRowsB, &rBeta, rResultC, &rNumRowsA );
            }
            else
            {
                //A*Bt
                if (rNumColumnsA!=rNumColumnsB)
                    throw NuTo::MathException("[FullMatrix::mult] number of columns of A must be equal to number of rows of B.");
                dgemm_ (&rTransA, &rTransB, &rNumRowsA, &rNumRowsB,  &rNumColumnsA, &rAlpha,
                        rA, &rNumRowsA, rB, &rNumRowsB, &rBeta, rResultC, &rNumRowsA );
            }
        }
        else
        {
            if (rTransB=='n')
            {
                //At*B
                if (rNumRowsA!=rNumRowsB)
                    throw NuTo::MathException("[FullMatrix::mult] number of columns of A must be equal to number of rows of B.");
                dgemm_ (&rTransA, &rTransB, &rNumColumnsA, &rNumColumnsB,  &rNumRowsA, &rAlpha,
                        rA, &rNumRowsA, rB, &rNumRowsB, &rBeta, rResultC, &rNumColumnsA );
            }
            else
            {
                //At*Bt
                if (rNumRowsA!=rNumColumnsB)
                    throw NuTo::MathException("[FullMatrix::mult] number of columns of A must be equal to number of rows of B.");
                dgemm_ (&rTransA, &rTransB, &rNumColumnsA, &rNumRowsB,  &rNumRowsA, &rAlpha,
                        rA, &rNumRowsA, rB, &rNumRowsB, &rBeta, rResultC, &rNumColumnsA );
            }
        }
    }
    catch (...)
    {
        throw NuTo::MathException("[FullMatrix::mult]Error calling the dgemm routine in the external library.");
    }
#else
    throw NuTo::MathException("[FullMatrix<double>::mult]NuTo wasn't compiled with LAPACK.");
#endif
}

void NuTo::MatrixOperations::print (const double *rA, unsigned int rNumRowsA, unsigned int rNumColumnsA)
{
    NuTo::MatrixOperations::print(rA,rNumRowsA,rNumColumnsA,10,15);
}

void NuTo::MatrixOperations::print (const double *rA, unsigned int rNumRows, unsigned int rNumColumns,
                                    unsigned int precision, unsigned int width)
{
    std::cout.precision(precision);
    for (unsigned int count=0; count<rNumRows; count++)
    {
        for (unsigned int count2=0; count2<rNumColumns; count2++)
        {
            std::ostringstream format_message;
            format_message << std::setprecision(precision) << std::setw(width) << std::showpoint << rA[count+count2*rNumRows];
            std::cout << std::setw(width) << format_message.str() <<" ";
        }
        std::cout<<std::endl;
    }
}
