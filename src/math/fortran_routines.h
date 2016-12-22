// $Id$

#pragma once

extern "C"
{
    int dgemm_(char *transa, const char *transb, int *m, int *n, int *k,
    double *alpha, const double *a, int *lda, const double *b, int *ldb,
    double *beta, double *c, int *ldc);
}
