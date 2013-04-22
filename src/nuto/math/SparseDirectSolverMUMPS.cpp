// $Id$
#ifdef SHOW_TIME
    #include <ctime>
#endif

# ifdef _OPENMP
#include <omp.h>
# endif

#ifdef HAVE_MUMPS
extern "C"
{
#include <dmumps_c.h>
}
#endif // HAVE_MUMPS

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseDirectSolver.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"

NuTo::SparseDirectSolverMUMPS::SparseDirectSolverMUMPS() : SparseDirectSolver()
{
#ifdef HAVE_MUMPS
    // set default solver parameters
    // this->orderingType = 2;          // set ordering to METIS
#else // HAVE_MUMPS
    throw NuTo::MathException("[SparseDirectSolverMUMPS::SparseDirectSolverMUMPS] MUMPS-solver was not found on your system (check cmake)");
#endif // HAVE_MUMPS
}

#ifdef HAVE_MUMPS
void NuTo::SparseDirectSolverMUMPS::Solve(const NuTo::SparseMatrixCSR<double>& rMatrix, const NuTo::FullVector<double,Eigen::Dynamic>& rRhs, NuTo::FullVector<double,Eigen::Dynamic>& rSolution)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    // check rMatrix
    if (rMatrix.HasZeroBasedIndexing())
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::solve] one based indexing of sparse matrix is required for this solver.");
    }
    int matrixDimension = rMatrix.GetNumRows();
    if (matrixDimension != rMatrix.GetNumColumns())
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::solve] matrix is not square.");
    }
    const std::vector<int>& matrixRowIndex = rMatrix.GetRowIndex();
    // extract rows from rowIndex
    std::vector<int> matrixRows(rMatrix.GetNumEntries());
    int entryCount = 0;
    for (int rowCount = 0; rowCount < matrixDimension; rowCount++)
    {
        for (int indexCount = matrixRowIndex[rowCount]; indexCount < matrixRowIndex[rowCount + 1]; indexCount++)
        {
            assert(entryCount < rMatrix.GetNumEntries());
            matrixRows[entryCount] = rowCount + 1;
            entryCount++;
        }
    }
    const std::vector<int>& matrixColumns = rMatrix.GetColumns();
    const std::vector<double>& matrixValues = rMatrix.GetValues();

    // check right hand side
    if (matrixDimension != rRhs.GetNumRows())
    {
        throw NuTo::MathException("[SparseDirectSolverMKLPardiso::solve] invalid dimension of right hand side vector.");
    }
    int rhsNumColumns = rRhs.GetNumColumns();

    // prepare rSolution rMatrix (copy rMatrix of right hand side vectors)
    rSolution = rRhs;
    const double *solutionValues = rSolution.data();

    // initialize solver data
    DMUMPS_STRUC_C solver;
    // set MPI communicator (also required in sequential version)
    solver.comm_fortran = -987654;
    // host is involved in factorization and rSolution phase
    solver.par = 1;
    // set rMatrix type
    if (rMatrix.IsSymmetric())
    {
        if (rMatrix.IsPositiveDefinite())
        {
            solver.sym = 1;
        }
        else
        {
            solver.sym = 2;
        }
    }
    else
    {
        solver.sym = 0;
    }
    // initialize solver
    solver.job = -1;
    dmumps_c(&solver);

    // define the problem
    solver.n   = matrixDimension;                       // dimension
    solver.nz  = rMatrix.GetNumEntries();                // number of nonzero entries
    solver.irn = &matrixRows[0];                        // rows
    solver.jcn = const_cast<int*>(&matrixColumns[0]);   // columns
    solver.a   = const_cast<double*>(&matrixValues[0]); // values
    solver.rhs = 0;                                     // right hand side vector is set befor solution
    // define solver specific parameters
    solver.icntl[0] = 0; // output stream for error messages
    solver.icntl[1] = 0; // output stream for diagnostic printing, statistics, and warning messages
    solver.icntl[2] = 0; // output stream for global information
    solver.icntl[3] = 0; // level of printing for error, warning, and diagnostic messages (0..4)
    solver.icntl[5] = 7; //  control an option for permuting and/or scaling the rMatrix (7 - automaticchoice)
    solver.icntl[6] = 7; // determines the pivot order to be used for the factorization (7 -  automatic choice)
    // set scaling strategy
    if (solver.sym == 0)
    {
        solver.icntl[7] = 7; // simultaneous row and colum iterative scaling
    }
    else
    {
        solver.icntl[7] = 1; // diagonal scaling
    }

    // analysis phase
    solver.job = 1;
    dmumps_c(&solver);
    if (solver.info[0] < 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::solve] Analysis and reordering phase: " + this->GetErrorString(solver.info[0]) + ".");
    }

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout<<"[NuTo::SparseDirectSolverMUMPS::Solve] reordering and symbolic factorization " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
    start = end;
#endif

    // Numerical factorization.
    solver.job = 2;
    dmumps_c(&solver);
    if (solver.info[0] < 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::solve] Numerical factorization phase: " + this->GetErrorString(solver.info[0]) + ".");
    }

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout<<"[NuTo::SparseDirectSolverMUMPS::Solve] numerical factorization " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
    start = end;
#endif

    // solution
    solver.job = 3;
    for (int rhsCount = 0; rhsCount < rhsNumColumns; rhsCount++)
    {
        solver.rhs = const_cast<double*>(&solutionValues[rhsCount*matrixDimension]);
        dmumps_c(&solver);
        if (solver.info[0] < 0)
        {
            throw NuTo::MathException("[SparseDirectSolverMUMPS::solve] Solution phase: " + this->GetErrorString(solver.info[0]) + ".");
        }
    }

    // Termination and release of memory
    solver.job = -2;
    dmumps_c(&solver);
    if (solver.info[0] < 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::solve] Termination phase: " + this->GetErrorString(solver.info[0]) + ".");
    }
#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout<<"[NuTo::SparseDirectSolverMUMPS::Solve] Solution " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

void NuTo::SparseDirectSolverMUMPS::SchurComplement(const NuTo::SparseMatrixCSR<double>& rMatrix, NuTo::FullVector<int, Eigen::Dynamic> rSchurIndices, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSchurComplement)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    // check rMatrix
    if (rMatrix.HasZeroBasedIndexing())
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::SchurComplement] one based indexing of sparse matrix is required for this solver.");
    }
    int matrixDimension = rMatrix.GetNumRows();
    if (matrixDimension != rMatrix.GetNumColumns())
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::SchurComplement] matrix is not square.");
    }
    const std::vector<int>& matrixRowIndex = rMatrix.GetRowIndex();
    // extract rows from rowIndex
    std::vector<int> matrixRows(rMatrix.GetNumEntries());
    int entryCount = 0;
    for (int rowCount = 0; rowCount < matrixDimension; rowCount++)
    {
        for (int indexCount = matrixRowIndex[rowCount]; indexCount < matrixRowIndex[rowCount + 1]; indexCount++)
        {
            assert(entryCount < rMatrix.GetNumEntries());
            matrixRows[entryCount] = rowCount + 1;
            entryCount++;
        }
    }
    const std::vector<int>& matrixColumns = rMatrix.GetColumns();
    const std::vector<double>& matrixValues = rMatrix.GetValues();

    // check Schur Indices
    if (matrixDimension < rSchurIndices.GetNumRows() || rSchurIndices.GetNumColumns()!=1)
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::SchurComplement] invalid dimension of schur indices vector.");
    }

    // add indices by one, since the external interface always counts from zero and mumps requires one based indexing
    rSchurIndices+=1;

    // initialize solver data
    DMUMPS_STRUC_C solver;
    // set MPI communicator (also required in sequential version)
    solver.comm_fortran = -987654;
    // host is involved in factorization and rSolution phase
    solver.par = 1;
    // set rMatrix type
    if (rMatrix.IsSymmetric())
    {
        if (rMatrix.IsPositiveDefinite())
        {
            solver.sym = 1;
        }
        else
        {
            solver.sym = 2;
        }
    }
    else
    {
        solver.sym = 0;
    }
    // initialize solver
    solver.job = -1;
    dmumps_c(&solver);

    //resize result matrix
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSchurComplementTranspose(rSchurIndices.GetNumRows(),rSchurIndices.GetNumRows());

    // define the problem
    solver.n   = matrixDimension;                       // dimension
    solver.nz  = rMatrix.GetNumEntries();               // number of nonzero entries
    solver.irn = &matrixRows[0];                        // rows
    solver.jcn = const_cast<int*>(&matrixColumns[0]);   // columns
    solver.a   = const_cast<double*>(&matrixValues[0]); // values
    solver.rhs = 0;
    solver.size_schur = rSchurIndices.GetNumRows();
    solver.listvar_schur = const_cast<int*>(rSchurIndices.data()); //variables of the schur matrix (the 1..size_schur)
    solver.schur = rSchurComplementTranspose.data();

    // define solver specific parameters
    solver.icntl[0] = 0; // output stream for error messages
    solver.icntl[1] = 0; // output stream for diagnostic printing, statistics, and warning messages
    solver.icntl[2] = 0; // output stream for global information
    solver.icntl[3] = 0; // level of printing for error, warning, and diagnostic messages (0..4)
    solver.icntl[5] = 0; // control an option for permuting and/or scaling the rMatrix (0 - no scaling due to Schur)
    solver.icntl[6] = 7; // determines the pivot order to be used for the factorization (7 -  automatic choice)
    solver.icntl[7] = 77; // set scaling strategy to automatic
    //solver.icntl[13] = 1000; // additional fill in (1000%, standard is 20%) this might indicate that you have not removed the almost zero entries of your matrix
    solver.icntl[18] = 1; // centralized Schur complement

    // analysis phase
    solver.job = 1;
    dmumps_c(&solver);
    if (solver.info[0] < 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::SchurComplement] Analysis and reordering phase: " + this->GetErrorString(solver.info[0]) + ".");
    }

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout<<"[NuTo::SparseDirectSolverMUMPS::SchurComplement] reordering and symbolic factorization " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
    start = end;
#endif

    // Numerical factorization.
    solver.job = 2;
    dmumps_c(&solver);
    if (solver.info[0] < 0)
    {
        std::cout << "solver info " << solver.info[0] << "\n";
    	throw NuTo::MathException("[SparseDirectSolverMUMPS::SchurComplement] Numerical factorization phase: " + this->GetErrorString(solver.info[0]) + ".");
    }

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout<<"[NuTo::SparseDirectSolverMUMPS::SchurComplement] numerical factorization " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
    start = end;
#endif

    // Termination and release of memory
    solver.job = -2;
    dmumps_c(&solver);
    if (solver.info[0] < 0)
    {
        throw NuTo::MathException("[SparseDirectSolverMUMPS::SchurComplement] Termination phase: " + this->GetErrorString(solver.info[0]) + ".");
    }

    rSchurComplement = rSchurComplementTranspose.transpose();

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout<<"[NuTo::SparseDirectSolverMUMPS::SchurComplement] Solution " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


std::string NuTo::SparseDirectSolverMUMPS::GetErrorString(int error) const
{
    assert(error < 0);

    switch (error)
    {
    case -9:
        return "main work array too small";
        break;
    case -10:
        return "numerically singular matrix";
        break;
    default:
        return "unknown error code";
    }
}
#else // HAVE_MUMPS
void NuTo::SparseDirectSolverMUMPS::Solve(const NuTo::SparseMatrixCSR<double>& rMatrix, const NuTo::FullVector<double, Eigen::Dynamic>& rRhs, NuTo::FullVector<double, Eigen::Dynamic>& rSolution)
{
	throw NuTo::MathException("[SparseDirectSolverMUMPS::Solve] MUMPS-solver was not found on your system (check cmake)");
}

void NuTo::SparseDirectSolverMUMPS::SchurComplement(const NuTo::SparseMatrixCSR<double>& rMatrix, const NuTo::FullVector<double, Eigen::Dynamic>& rRhs, NuTo::FullVector<double, Eigen::Dynamic>& rSolution)
{
    throw NuTo::MathException("[SparseDirectSolverMUMPS::SchurComplement] MUMPS-solver was not found on your system (check cmake)");
}
#endif // HAVE_MUMPS
