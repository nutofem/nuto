#pragma once

#ifdef HAVE_MUMPS
extern "C" {
#include <dmumps_c.h>
}
#endif // HAVE_MUMPS


#include "math/SparseDirectSolver.h"
#include <Eigen/Core>

namespace NuTo
{
// forward declarations
template <class T>
class SparseMatrixCSR;
class EigenSolverArpack;

//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... interface for the MUMPS sparse direct solver
class SparseDirectSolverMUMPS : public SparseDirectSolver
{
    friend class NuTo::EigenSolverArpack;

public:
    //! @brief ... default constructor
    SparseDirectSolverMUMPS();

    //! @brief ... solve system of equations: rMatrix * rSolution = rRhs
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    //! @param rRhs ... matrix storing the right-hand-side vectors (input)
    //! @param rSolution ... matrix storing the corresponding solution vectors (output)
    void Solve(const NuTo::SparseMatrixCSR<double>& rMatrix, const Eigen::VectorXd& rRhs, Eigen::VectorXd& rSolution);


    //! @brief ... calculates the Schurcomplement
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    //! @param rSchurIndices ... vector storing the indices of the global matrix to be condensed to (zero based
    //! indexing)
    //! @param rSchurComplement ... Schur complement
    void SchurComplement(const NuTo::SparseMatrixCSR<double>& rMatrix, Eigen::VectorXi rSchurIndices,
                         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSchurComplement);

    //! @brief ... prepare the solver, and perform all the steps up to the factorization of the matrix
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    void Factorization(const NuTo::SparseMatrixCSR<double>& rMatrix);

    //! @brief ... use the factorized matrix for the final solution phase
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    //! @param rRhs ... matrix storing the right-hand-side vectors (input)
    //! @param rSolution ... matrix storing the corresponding solution vectors (output)
    void Solution(const Eigen::VectorXd& rRhs, Eigen::VectorXd& rSolution);

    //! @brief ... Termination and release of memory
    void CleanUp();

protected:
#ifdef HAVE_MUMPS

    DMUMPS_STRUC_C mSolver;

#endif // HAVE_MUMPS

    //! @brief ... generate an error message from the error code
    //! @param error ... error code
    //! @return error message as std::string
    std::string GetErrorString(int error) const;
};
}
