// $Id$

#include "nuto/math/fortran_routines.h"


void mult(double *rA, int rNumRowsA, int rNumColumnsA,
          double *rB, int rNumRowsB, int rNumColumnsB,
          char rTransA, char rTransB, double rAlpha,
          double *rResultC, double rBeta )
{
#ifdef HAVE_LAPACK

    try
    {
        if (rTransA=='n')
        {
            if (rTransB=='n')
            {
                if (rNumColumnsA!=rNumRowsB)
                    throw MathException("[FullMatrix::mult] number of columns of A must be equal to number of rows of B.");
                dgemm_ ("n", "n", &rNumRowsA, &rNumColumnsB,  &rNumColumnsA, &rAlpha,
                        rA, &rNumRowsA, rB, &rNumColumnsA, &rBeta, rResultC, &rNumRowsA );
            }
            else
            {
                throw MathException("[FullMatrix::mult]Not yet implemented.");
            }
        }
    }
    catch (...)
    {
        throw MathException("[FullMatrix::mult]Error calling the dgemm routine in the external library.");
    }
#else
    throw MathException("[FullMatrix<double>::mult]NuTo wasn't compiled with LAPACK.");
#endif
}
