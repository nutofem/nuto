#pragma once

#include "math/SparseDirectSolver.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{
// forward declarations
template<class T> class SparseMatrixCSR;

//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... interface for the MKL-DSS sparse direct solver
class SparseDirectSolverMKLDSS  : public SparseDirectSolver
{
public:
    //! @brief ... default constructor
    SparseDirectSolverMKLDSS();

    //! @brief ... solve the system of equations
    //! @param rMatrix ... sparse coefficient matrix stored in the CSR format
    //! @param rRhs ... matrix of right-hand-side vectors (full storage)
    //! @param rSolution ... matrix of solution vectors (full storage)
    void Solve(const SparseMatrixCSR<double>& rMatrix, const Eigen::VectorXd& rRhs, Eigen::VectorXd& rSolution);

#ifdef HAVE_MKL_DSS

    //! @brief ... enable refinement in the solution stage
    inline void enableRefinement()
    {
        this->mRefinement = true;
    }

    //! @brief ... disable refinement in the solution stage
    inline void DisableRefinement()
    {
        this->mRefinement = false;
    }
protected:
    //! @brief ... determines wether a refinement is enabled or disabled within the solution strategy
    bool mRefinement;
#endif // HAVE_MKL_DSS
};
}
