// $Id$

#ifndef NUTOMATH_H
#define NUTOMATH_H

namespace NuTo
{
typedef enum
{
    BINARY=0,
    XML,
    TEXT
} serialization_attributes;

class MatrixOperations
{
public:
#ifndef SWIG
    static void mult(const double *rA, int rNumRowsA, int rNumColumnsA,
                     const double *rB, int rNumRowsB, int rNumColumnsB,
                     char rTransA, char rTransB, double rAlpha,
                     double *rResultC, double rBeta);

    static void sub (const double *rA, unsigned int rNumRowsA, unsigned int rNumColumnsA,
                     const double *rB, double *rC)
    {
        for (unsigned int count=0; count<rNumRowsA*rNumColumnsA; count++,rC++,rA++,rB++)
        {
            *rC=*rA-(*rB);
        }
    }

    static void print (const double *rA, unsigned int rNumRowsA, unsigned int rNumColumnsA,
                       unsigned int precision, unsigned int width);
    static void print (const double *rA, unsigned int rNumRowsA, unsigned int rNumColumnsA);
#endif
};
}
#endif // NUTOMATH_H
