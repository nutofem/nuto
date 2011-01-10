// $Id$

// header for lapack routines

#ifdef __cplusplus
extern "C" {
#endif

//! @brief DPOTRI computes the inverse of a real symmetric positive definite matrix A using the Cholesky factorization A = U**T*U or A = L*L**T computed by DPOTRF.
//! @param UPLO ... (input) CHARACTER*1, 'U':  Upper triangle of A is stored, 'L':  Lower triangle of A is stored.
//! @param N ... (input) INTEGER, The order of the matrix A.  N >= 0.
//! @param A ... (input/output) DOUBLE PRECISION array, dimension (LDA,N), On entry, the triangular factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T, as computed by DPOTRF. On exit, the upper or lower triangle of the (symmetric) inverse of A, overwriting the input factor U or L.
//! @param LDA ... (input) INTEGER, The leading dimension of the array A.  LDA >= max(1,N).
//! @param INFO    (output) INTEGER, = 0:  successful exit, < 0:  if INFO = -i, the i-th argument had an illegal value, > 0:  if INFO = i, the (i,i) element of the factor U or L is zero, and the inverse could not be computed.
void dpotri_(const char* UPLO, const int* N, double* A, const int* LDA, int* INFO );

//! @brief DPOTRF computes the Cholesky factorization of a real symmetric positive definite matrix A. The factorization has the form A = U**T * U,  if UPLO = 'U', or A = L  * L**T,  if UPLO = 'L', *  where U is an upper triangular matrix and L is lower triangular. This is the block version of the algorithm, calling Level 3 BLAS.
//! @param UPLO ... (input) CHARACTER*1, 'U':  Upper triangle of A is stored, 'L':  Lower triangle of A is stored.
//! @param N ... (input) INTEGER, The order of the matrix A.  N >= 0.
//! @param A ... (input/output) DOUBLE PRECISION array, dimension (LDA,N), On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the upper triangular part of the matrix A, and the strictly lower triangular part of A is not referenced.  If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A, and the strictly upper triangular part of A is not referenced. On exit, if INFO = 0, the factor U or L from the Cholesky factorization A = U**T*U or A = L*L**T.
//! @param LDA ... (input) INTEGER, The leading dimension of the array A.  LDA >= max(1,N).
//! @param INFO ... (output) INTEGER, = 0:  successful exit, < 0:  if INFO = -i, the i-th argument had an illegal value, > 0:  if INFO = i, the leading minor of order i is not positive definite, and the factorization could not be completed.
void dpotrf_(const char* UPLO, const int* N, double* A, const int* LDA, int* INFO);

//! @brief DPOTRS solves a system of linear equations A*X = B with a symmetric positive definite matrix A using the Cholesky factorization A = U**T*U or A = L*L**T computed by DPOTRF.
//! @param UPLO ... (input) CHARACTER*1, 'U':  Upper triangle of A is stored, 'L':  Lower triangle of A is stored.
//! @param N ... (input) INTEGER, The order of the matrix A.  N >= 0.
//! @param NRHS ... (input) INTEGER, The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
//! @param A ... (input) DOUBLE PRECISION array, dimension (LDA,N), The triangular factor U or L from the Cholesky factorization, A = U**T*U or A = L*L**T, as computed by DPOTRF.
//! @param LDA ... (input) INTEGER, The leading dimension of the array A.  LDA >= max(1,N).
//! @param B ... (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS), On entry, the right hand side matrix B. On exit, the solution matrix X.
//! @param LDB ... (input) INTEGER, The leading dimension of the array B.  LDB >= max(1,N).
//! @param INFO ... (output) INTEGER, = 0:  successful exit, < 0:  if INFO = -i, the i-th argument had an illegal value
void dpotrs_(const char* UPLO, const int* N, const int* NRHS, const double* A, const int* LDA, double* B, const int* LDB, int* INFO);

#ifdef __cplusplus
}
#endif
