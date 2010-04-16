#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H
namespace NuTo
{
//! @author Joerg F. Unger, ISM
//! @date July 2009
//! @brief ... class for full matrices derived from the abstract base class Matrix
template <class T>
class MatrixOperations
{
public:
    static void mult(double *rA, int rNumRowsA, int rNumColumnsA,
                     double *rB, int rNumRowsB, int rNumColumnsB,
                     char rTransA, char rTransB, double rAlpha,
                     double *rResultC, double rBeta);

    static void sub (double *rA, int rNumRowsA, int rNumColumnsA,
                     double *rB, double *rC)
    {
        for (unsigned int count=0; count<rNumRowsA*rNumColumnsA; count++)
        {
            *rC=*rA-(*rB);
            rC++;
            rA++;
            rB++;
        }
    }
};
}
#endif MATRIX_OPERATIONS_H
