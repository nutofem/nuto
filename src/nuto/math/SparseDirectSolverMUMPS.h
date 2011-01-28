// $Id$
#ifndef SPARSE_DIRECT_SOLVER_MUMPS_H
#define SPARSE_DIRECT_SOLVER_MUMPS_H

#include "nuto/math/SparseDirectSolver.h"

namespace NuTo
{
// forward declarations
template<class T> class SparseMatrixCSR;
template<class T> class FullMatrix;

//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... interface for the MUMPS sparse direct solver
class SparseDirectSolverMUMPS : public SparseDirectSolver
{
public:
    //! @brief ... default constructor
    SparseDirectSolverMUMPS();

    //! @brief ... print information about the class attributes
    virtual void Info()const
    {
    }

    std::string GetTypeId()const
    {
        return std::string("SparseDirectSolverMUMPS");
    }

#ifdef HAVE_MUMPS
    //! @brief ... solve system of equations: rMatrix * rSolution = rRhs
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    //! @param rRhs ... matrix storing the right-hand-side vectors (input)
    //! @param rSolution ... matrix storing the corresponding solution vectors (output)

    void Solve(const SparseMatrixCSR<double>& rMatrix, const FullMatrix<double>& rRhs, FullMatrix<double>& rSolution);

    //! @brief ... calculates the Schurcomplement
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    //! @param rSchurIndices ... vector/matrix storing the indices of the global matrix to be condensed to (zero based indexing)
    //! @param rSchurComplement ... Schur complement
    void SchurComplement(const NuTo::SparseMatrixCSR<double>& rMatrix, NuTo::FullMatrix<int> rSchurIndices, NuTo::FullMatrix<double>& rSchurComplement);

protected:
    //! @brief ... generate an error message from the error code
    //! @param error ... error code
    //! @return error message as std::string
    std::string GetErrorString(int error) const;

    //! @brief ... type of fill-in reducing odering of the coefficient matrix
    //int orderingType;
#else
    //! @brief ... solve system of equations: rMatrix * rSolution = rRhs
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    //! @param rRhs ... matrix storing the right-hand-side vectors (input)
    //! @param rSolution ... matrix storing the corresponding solution vectors (output)
    void Solve(const SparseMatrixCSR<double>& rMatrix, const FullMatrix<double>& rRhs, FullMatrix<double>& rSolution);

    void SchurComplement(const NuTo::SparseMatrixCSR<double>& rMatrix, const NuTo::FullMatrix<double>& rRhs, NuTo::FullMatrix<double>& rSolution);

#endif // HAVE_MUMPS
};
}
#endif // SPARSE_DIRECT_SOLVER_MUMPS_H
