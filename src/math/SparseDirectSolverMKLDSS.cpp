
#include <ctime>

#include "base/Exception.h"
#include "math/SparseMatrixCSR.h"
#include "math/SparseDirectSolver.h"
#include "math/SparseDirectSolverMKLDSS.h"

#ifdef HAVE_MKL_DSS
#include <mkl_dss.h>
#endif

NuTo::SparseDirectSolverMKLDSS::SparseDirectSolverMKLDSS()
    : SparseDirectSolver()
{
#ifdef HAVE_MKL_DSS
    this->mRefinement = true;
#else // HAVE_MKL_DSS
    throw NuTo::Exception("MKL DSS-solver was not found on your system (check cmake)");
#endif // HAVE_MKL_DSS
}

#ifdef HAVE_MKL_DSS
void NuTo::SparseDirectSolverMKLDSS::Solve(const NuTo::SparseMatrixCSR<double>& rMatrix, const Eigen::VectorXd& rRhs,
                                           Eigen::VectorXd& rSolution)
{
    // timing
    clock_t startTime = clock();

    // check rMatrix
    if (rMatrix.HasZeroBasedIndexing())
    {
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve] one based indexing of sparse matrix is required for this solver.");
    }
    int matrixDimension = rMatrix.GetNumRows();
    if (matrixDimension != rMatrix.GetNumColumns())
    {
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve] matrix must be symmetric.");
    }
    const std::vector<int>& matrixRowIndex = rMatrix.GetRowIndex();
    const std::vector<int>& matrixColumns = rMatrix.GetColumns();
    const std::vector<double>& matrixValues = rMatrix.GetValues();
    int matrixNumEntries = rMatrix.GetNumEntries();

    // check right hand side
    if (matrixDimension != rRhs.GetNumRows())
    {
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve] invalid dimension of right hand side vector.");
    }
    int rhsNumColumns = rRhs.GetNumColumns();
    const double* rhsValues = rRhs.GetEigenMatrix().data();

    // prepare solution matrix
    rSolution.Resize(matrixDimension, rhsNumColumns);
    const double* solutionValues = rSolution.GetEigenMatrix().data();

    // initialize solver
    _MKL_DSS_HANDLE_t handle;
    int options = MKL_DSS_DEFAULTS;
    int error = dss_create(handle, options);
    switch (error)
    {
    case MKL_DSS_SUCCESS:
        break;
    case MKL_DSS_INVALID_OPTION:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_create: invalid options.");
    case MKL_DSS_OUT_OF_MEMORY:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_create: out of memory.");
    default:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_create: unknown error code.");
    }

    // communicate rMatrix structure to solver
    int symmetry;
    if (rMatrix.IsSymmetric())
    {
        symmetry = MKL_DSS_SYMMETRIC;
    }
    else
    {
        symmetry = MKL_DSS_NON_SYMMETRIC;
    }
    error = dss_define_structure(handle, symmetry, &matrixRowIndex[0], matrixDimension, matrixDimension,
                                 &matrixColumns[0], matrixNumEntries);
    switch (error)
    {
    case MKL_DSS_SUCCESS:
        break;
    case MKL_DSS_STATE_ERR:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_define_structure: state error.");
    case MKL_DSS_INVALID_OPTION:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_define_structure: invalid options.");
    case MKL_DSS_COL_ERR:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_define_structure: column error.");
    case MKL_DSS_NOT_SQUARE:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_define_structure: matrix is not square.");
    case MKL_DSS_TOO_FEW_VALUES:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_define_structure: matrix has too few entries.");
    case MKL_DSS_TOO_MANY_VALUES:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_define_structure: matrix has too many entries.");
    default:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_define_structure: unknown error code.");
    }

    // reordering
    int reorder = MKL_DSS_AUTO_ORDER;
    error = dss_reorder(handle, reorder, NULL);
    switch (error)
    {
    case MKL_DSS_SUCCESS:
        break;
    case MKL_DSS_STATE_ERR:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_reorder: state error.");
    case MKL_DSS_INVALID_OPTION:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_reorder: invalid options.");
    case MKL_DSS_OUT_OF_MEMORY:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_reorder: out of memory.");
    default:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_reorder: unknown error code.");
    }
    if (this->mVerboseLevel > 0)
    {
        clock_t endTime = clock();
        std::cout << "[SparseDirectSolverMKLDSS::solve] Time for reordering and symbolic factorization: "
                  << static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
        startTime = endTime;
    }

    // factorization
    int type;
    if (rMatrix.IsPositiveDefinite())
    {
        type = MKL_DSS_POSITIVE_DEFINITE;
    }
    else
    {
        type = MKL_DSS_INDEFINITE;
    }
    error = dss_factor_real(handle, type, &matrixValues[0]);
    switch (error)
    {
    case MKL_DSS_SUCCESS:
        break;
    case MKL_DSS_STATE_ERR:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_factor_real: state error.");
    case MKL_DSS_INVALID_OPTION:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_factor_real: invalid options.");
    case MKL_DSS_OPTION_CONFLICT:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_factor_real: option conflict.");
    case MKL_DSS_OUT_OF_MEMORY:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_factor_real: out of memory.");
    case MKL_DSS_ZERO_PIVOT:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_factor_real: zero pivot.");
    default:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_factor_real: unknown error code.");
    }
    if (this->mVerboseLevel > 0)
    {
        clock_t endTime = clock();
        std::cout << "[SparseDirectSolverMKLDSS::solve] Time for numerical factorization: "
                  << static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
        startTime = endTime;
    }

    // solve
    int refinementOptions;
    if (this->mRefinement)
    {
        refinementOptions = MKL_DSS_REFINEMENT_ON;
    }
    else
    {
        refinementOptions = MKL_DSS_REFINEMENT_OFF;
    }
    error = dss_solve_real(handle, refinementOptions, rhsValues, rhsNumColumns, const_cast<double*>(solutionValues));
    switch (error)
    {
    case MKL_DSS_SUCCESS:
        break;
    case MKL_DSS_STATE_ERR:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_solve_real: state error.");
    case MKL_DSS_INVALID_OPTION:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_solve_real: invalid options.");
    case MKL_DSS_OUT_OF_MEMORY:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_solve_real: out of memory.");
    default:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_solve_real: unknown error code.");
    }
    if (this->mVerboseLevel > 0)
    {
        clock_t endTime = clock();
        std::cout << "[SparseDirectSolverMKLDSS::solve] Time for back substitution: "
                  << static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
        startTime = endTime;
    }

    if (this->mVerboseLevel > 1)
    {
        double statOut[4];
        // char statIn[] = "ReorderTime,FactorTime,SolveTime,Flops,Peakmem,Factormem,Solvemem";
        error = dss_statistics(handle, options, "Flops,Peakmem,Factormem,Solvemem", statOut);
        switch (error)
        {
        case MKL_DSS_SUCCESS:
            break;
        case MKL_DSS_STATISTICS_INVALID_MATRIX:
            throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: invalid matrix.");
        case MKL_DSS_STATISTICS_INVALID_STATE:
            throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: invalid state.");
        case MKL_DSS_STATISTICS_INVALID_STRING:
            throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: invalid string.");
        default:
            throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: unknown error code.");
        }
        std::cout << "[SparseDirectSolverMKLDSS::solve] Peak memory symbolic factorization: " << statOut[1] << " KBytes"
                  << std::endl;
        std::cout << "[SparseDirectSolverMKLDSS::solve] Permanent memory symbolic factorization: " << statOut[2]
                  << " KBytes" << std::endl;
        std::cout << "[SparseDirectSolverMKLDSS::solve] Memory numerical factorization and solution: " << statOut[3]
                  << " KBytes" << std::endl;
        if (this->mVerboseLevel > 2)
        {
            std::cout << "[SparseDirectSolverMKLDSS::solve] Number of floating point operations required for "
                         "factorization: "
                      << statOut[0] << " MFLOS" << std::endl;
            if (symmetry == MKL_DSS_SYMMETRIC && type == MKL_DSS_INDEFINITE)
            {
                error = dss_statistics(handle, options, "Inertia", statOut);
                switch (error)
                {
                case MKL_DSS_SUCCESS:
                    break;
                case MKL_DSS_STATISTICS_INVALID_MATRIX:
                    throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: invalid matrix.");
                case MKL_DSS_STATISTICS_INVALID_STATE:
                    throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: invalid state.");
                case MKL_DSS_STATISTICS_INVALID_STRING:
                    throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: invalid string.");
                default:
                    throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_statistics: unknown error code.");
                }
                std::cout << "[SparseDirectSolverMKLDSS::solve] Inertia: number of positive eigenvalues: " << statOut[0]
                          << std::endl;
                std::cout << "[SparseDirectSolverMKLDSS::solve] Inertia: number of negative eigenvalues: " << statOut[1]
                          << std::endl;
                std::cout << "[SparseDirectSolverMKLDSS::solve] Inertia: number of zero eigenvalues: " << statOut[2]
                          << std::endl;
            }
        }
    }

    error = dss_delete(handle, options);
    switch (error)
    {
    case MKL_DSS_SUCCESS:
        break;
    case MKL_DSS_INVALID_OPTION:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_delete: invalid options .");
    case MKL_DSS_OUT_OF_MEMORY:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_delete: out of memory.");
    default:
        throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]dss_delete: unknown error code.");
    }
}
#else // HAVE_MKL_DSS
void NuTo::SparseDirectSolverMKLDSS::Solve(const NuTo::SparseMatrixCSR<double>& rMatrix, const Eigen::VectorXd& rRhs,
                                           Eigen::VectorXd& rSolution)
{
	throw NuTo::Exception("[SparseDirectSolverMKLDSS::solve]MKLDSS not implemented.");
}
#endif // HAVE_MKL_DSS
