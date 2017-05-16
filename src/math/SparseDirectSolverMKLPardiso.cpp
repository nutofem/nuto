// $Id$

#include <ctime>

#ifdef HAVE_MKL_PARDISO
#include <mkl_pardiso.h>
#include <mkl_service.h>
#endif // HAVE_MKL_PARDISO

#include "base/Timer.h"
#include "math/MathException.h"
#include "math/SparseMatrixCSR.h"
#include "math/SparseDirectSolver.h"
#include "math/SparseDirectSolverMKLPardiso.h"


NuTo::SparseDirectSolverMKLPardiso::SparseDirectSolverMKLPardiso() : SparseDirectSolver()
{
#ifdef HAVE_MKL_PARDISO
    // set default solver parameters
    this->mOrderingType = 2;          // set ordering to METIS
    this->mNumRefinementSteps  = 0;   // maximum number of iterative refinement steps
    this->mPivotingPerturbation = -1; // set pivoting perturbation to default values
    this->mScaling = -1;              // enable nonsymmetric permutation and mScaling MPS only for unsymmetric matrices (MKL default)
    this->mWeightedMatching = -1;     // enable maximum weighted matching algorithm only for unsymmetric matrices (MKL default)
#else // HAVE_MKL_PARDISO
    throw NuTo::MathException("MKL Pardiso-solver was not found on your system (check cmake)");
#endif // HAVE_MKL_PARDISO
}

#ifdef HAVE_MKL_PARDISO
void NuTo::SparseDirectSolverMKLPardiso::Solve(const NuTo::SparseMatrixCSR<double>& rMatrix, const Eigen::VectorXd& rRhs, Eigen::VectorXd& rSolution)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__ + " Reordering and symbolic factorization");
    
    mkl_set_num_threads(2);
    std::cout << "number of threads: " << mkl_get_max_threads() << std::endl;

    // check rMatrix
    if (rMatrix.HasZeroBasedIndexing())
    {
        throw NuTo::MathException(__PRETTY_FUNCTION__, "one based indexing of sparse rMatrix is required for this solver.");
    }
    int matrixDimension = rMatrix.GetNumRows();
    if (matrixDimension != rMatrix.GetNumColumns())
    {
        throw NuTo::MathException(__PRETTY_FUNCTION__, "matrix must be symmetric.");
    }
    const std::vector<int>& matrixRowIndex = rMatrix.GetRowIndex();
    const std::vector<int>& matrixColumns = rMatrix.GetColumns();
    const std::vector<double>& matrixValues = rMatrix.GetValues();
    int matrixType;
    if (rMatrix.IsSymmetric())
    {
        if (rMatrix.IsPositiveDefinite())
        {
            matrixType = 2;
        }
        else
        {
            matrixType = -2;
        }
    }
    else
    {
        matrixType = 11;
    }

    // check right hand side
    if (matrixDimension != rRhs.rows())
    {
        throw NuTo::MathException(__PRETTY_FUNCTION__, "invalid dimension of right hand side vector.");
    }
    int rhsNumColumns = rRhs.cols();
    const double *rhsValues = rRhs.data();

    // prepare solution matrix
    rSolution.Resize(matrixDimension,rhsNumColumns);
    const double *solutionValues = rSolution.data();

    // initialize solver data
    int parameters[64];
    for (unsigned int parameterCount = 0; parameterCount < 64; parameterCount++)
    {
        parameters[parameterCount] = 0;
    }
    // use non-default solver parameters
    parameters[0] = 1;
    // set ordering to METIS
    parameters[1] = this->mOrderingType;
    // set number of threads
    parameters[2] = 2;
    // set maximum number of iterative refinement steps
    parameters[7] = this->mNumRefinementSteps;
    // set perturbation of pivot elements
    if (this->mPivotingPerturbation < 0)
    {
        // default values
        if (matrixType == 11)
        {
            parameters[9] = 13;
        }
        else
        {
            parameters[9] = 8;
        }
    }
    else
    {
        parameters[9] = this->mPivotingPerturbation;
    }
    // enable/disable nonsymmetric permutation and mScaling MPS
    if (this->mScaling < 0)
    {
        // use default
        if (matrixType == 11)
        {
            parameters[10] = 1;
        }
        else
        {
            parameters[10] = 0;
        }
    }
    else
    {
        parameters[10] = this->mScaling;
    }
    // enable/disable maximum weighted matching algorithm
    if (this->mWeightedMatching < 0)
    {
        // use default
        if (matrixType == 11)
        {
            parameters[12] = 1;
        }
        else
        {
            parameters[12] = 0;
        }
    }
    else
    {
        parameters[12] = this->mWeightedMatching;
    }
    // enable additional output
    if (this->mVerboseLevel > 2)
    {
        // determine numbers of non-zero elements on the factors
        parameters[17] = -1;
        // determine number of MFlops (1.0E6) that are necessary to factor the rMatrix.
        // This will increase the reordering time.
        parameters[18] = -1;
    }
    void* pt[64];
    for (unsigned int count = 0; count < 64; count++)
    {
        pt[count] = 0;
    }
    int maxfct(1);  // Maximum number of numerical factorizations.
    int mnum(1);    // Which factorization to use.
    int msglvl(0);  // Print statistical information in file
    int error(0);   // Initialize error flag
    double ddum(0); // Double dummy
    int idum(0);    // Integer dummy

    // Reordering and Symbolic Factorization.
    // This step also allocates all memory that is necessary for the factorization.
    int phase = 11;
    PARDISO (pt,
             &maxfct,
             &mnum,
             &matrixType,
             &phase,
             &matrixDimension,
             const_cast<double*>(&matrixValues[0]),
             const_cast<int*>(&matrixRowIndex[0]),
             const_cast<int*>(&matrixColumns[0]),
             &idum,
             &rhsNumColumns,
             parameters,
             &msglvl,
             &ddum,
             &ddum,
             &error);
    if (error != 0)
    {
        throw NuTo::MathException(__PRETTY_FUNCTION__, "Analysis and reordering phase: " + this->GetErrorString(error) + ".");
    }
    
    timer.Restart(__PRETTY_FUNCTION__ + " Numerical factorization");

    phase = 22;
    PARDISO (pt,
             &maxfct,
             &mnum,
             &matrixType,
             &phase,
             &matrixDimension,
             const_cast<double*>(&matrixValues[0]),
             const_cast<int*>(&matrixRowIndex[0]),
             const_cast<int*>(&matrixColumns[0]),
             &idum,
             &rhsNumColumns,
             parameters,
             &msglvl,
             &ddum,
             &ddum,
             &error);
    if (error != 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMKLPardiso::solve] Numerical factorization phase: " + this->GetErrorString(error) + ".");
    }
    
    timer.Reset(__PRETTY_FUNCTION__ + "Back substitution and iterative refinement.");

    phase = 33;
    PARDISO (pt,
             &maxfct,
             &mnum,
             &matrixType,
             &phase,
             &matrixDimension,
             const_cast<double*>(&matrixValues[0]),
             const_cast<int*>(&matrixRowIndex[0]),
             const_cast<int*>(&matrixColumns[0]),
             &idum,
             &rhsNumColumns,
             parameters,
             &msglvl,
             const_cast<double*>(rhsValues),
             const_cast<double*>(solutionValues),
             &error);
    if (error != 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMKLPardiso::solve] Back substitution and iterative refinement phase: " + this->GetErrorString(error) + ".");
    }

	if (this->mVerboseLevel > 1)
	{
		std::cout << "[SparseDirectSolverMKLPardiso::solve] Peak memory symbolic factorization: "
				  << parameters[14] << " KBytes"
				  << std::endl;
		std::cout << "[SparseDirectSolverMKLPardiso::solve] Permanent memory symbolic factorization: "
				  << parameters[15] << " KBytes"
				  << std::endl;
		std::cout << "[SparseDirectSolverMKLPardiso::solve] Memory numerical factorization and solution: "
				  << parameters[16] << " KBytes"
				  << std::endl;
		if (this->mVerboseLevel > 2)
		{
			std::cout << "[SparseDirectSolverMKLPardiso::solve] Number of floating point operations required for factorization: "
					  << parameters[18] << " MFLOS"
					  << std::endl;
			if (matrixType == -2)
			{
				std::cout << "[SparseDirectSolverMKLPardiso::solve] Inertia: number of positive eigenvalues: "
						  << parameters[21]
						  << std::endl;
				std::cout << "[SparseDirectSolverMKLPardiso::solve] Inertia: number of negative eigenvalues: "
						  << parameters[22]
						  << std::endl;
				std::cout << "[SparseDirectSolverMKLPardiso::solve] Inertia: number of zero eigenvalues: "
						  << matrixDimension - parameters[21] - parameters[22]
						  << std::endl;
			}
			std::cout << "[SparseDirectSolverMKLPardiso::solve] Number of nonzeros in factors: "
					  << parameters[17]
					  << std::endl;
			std::cout << "[SparseDirectSolverMKLPardiso::solve] Number of performed iterative refinement steps: "
					  << parameters[6]
					  << std::endl;
			if (matrixType != 2)
			{
				std::cout << "[SparseDirectSolverMKLPardiso::solve] Number of perturbed pivots: "
						  << parameters[13]
						  << std::endl;
			}
		}
    }

    // Termination and release of memory
    phase = -1;
    PARDISO (pt,
             &maxfct,
             &mnum,
             &matrixType,
             &phase,
             &matrixDimension,
             &ddum,
             const_cast<int*>(&matrixRowIndex[0]),
             const_cast<int*>(&matrixColumns[0]),
             &idum,
             &rhsNumColumns,
             parameters,
             &msglvl,
             &ddum,
             &ddum,
             &error);
    if (error != 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMKLPardiso::solve] Termination phase: " + this->GetErrorString(error) + ".");
    }
}

std::string NuTo::SparseDirectSolverMKLPardiso::GetErrorString(int error) const
{
    assert(error != 0);

    switch (error)
    {
    case -1:
        return "input inconsistent";
    case -2:
        return "not enough memory";
    case -3:
        return "reordering problem";
    case -4:
        return "zero pivot, numerical factorization or iterative refinement problem";
    case -5:
        return "unclassified (internal) error";
    case -6:
        return "preordering failed";
    case -7:
        return "diagonal rMatrix problem";
    case -8:
        return "32-bit integer overflow problem";
    case -9:
        return "not enough memory for OOC";
    case -10:
        return "problems with opening OOC temporary files";
    case -11:
        return "read/write problems with the OOC data file";
    default:
        return "unknown error code";
    }
}
#endif // HAVE_MKL_PARDISO
