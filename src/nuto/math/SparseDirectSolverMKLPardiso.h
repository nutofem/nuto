// $Id$
#ifndef SPARSE_DIRECT_SOLVER_MKLPARDISO_H
#define SPARSE_DIRECT_SOLVER_MKLPARDISO_H

#include "nuto/math/SparseDirectSolver.h"

namespace NuTo
{
// forward declarations
template<class T> class SparseMatrixCSR;
template<class T> class FullMatrix;

//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... interface for the MKL-PARDISO sparse direct solver
class SparseDirectSolverMKLPardiso : public SparseDirectSolver
{
public:
    //! @brief ... default constructor
    SparseDirectSolverMKLPardiso();

    //! @brief ... print information about the class attributes
    virtual void Info()
    {
    }
#ifdef HAVE_MKL_PARDISO
    //! @brief ... solve system of equations: rMatrix * rSolution = rRhs
    //! @param rMatrix ... sparse coefficient matrix, stored in compressed CSR format (input)
    //! @param rRhs ... matrix storing the right-hand-side vectors (input)
    //! @param rSolution ... matrix storing the corresponding solution vectors (output)
    void Solve(const SparseMatrixCSR<double>& rMatrix, const FullMatrix<double>& rRhs, FullMatrix<double>& rSolution);

    //! @brief ... use the nested dissection alogrithm from the METIS-package for for the fill-in reducing odering of the coefficient matrix
    //! @sa mOrderingType
    inline void SetOrderingMETIS()
    {
        this->mOrderingType = 2;  // set ordering to METIS
    }

    //! @brief ... use the minimum degree alogrithm for for the fill-in reducing odering of the coefficient matrix
    //! @sa mOrderingType
    inline void SetOrderingMinimumDegree()
    {
        this->mOrderingType = 0;  // set ordering to METIS
    }

    //! @brief ... set the maximum number of iterative refinement steps aplied in the solution phase
    //! @param rNumRefinementSteps_ ... the number of maximum iterative refinement steps
    //! @sa mNumRefinementSteps
    inline void SetNumRefinementSteps(unsigned int rNumRefinementSteps_)
    {
        this->mNumRefinementSteps = rNumRefinementSteps_;
    }

    //! @brief ... set the perturbation for small pivots. Small pivots will be perturbed with epsilon = 10^(-orderOfMagnitude).
    //! @param rOrderOfMagnitude ... absolute value of the order of magnitude of the perturbation for small pivots
    //! @sa mPivotingPerturbation
    inline void SetPivotingPerturbation(unsigned int rOrderOfMagnitude)
    {
        this->mPivotingPerturbation = rOrderOfMagnitude;
    }

    //! @brief ... enable scaling of the coefficient matrix matrix
    //! @sa mScaling
    inline void EnableScaling()
    {
        this->mScaling = 1;
    }

    //! @brief ... disable scaling of the coefficient matrix matrix
    //! @sa mScaling
    inline void DisableScaling()
    {
        this->mScaling = 0;
    }

    //! @brief ... enable maximum weighted matching algorithm
    //! @sa mWeightedMatching
    inline void EnableWeightedMatching()
    {
        this->mWeightedMatching = 1;
    }

    //! @brief ... disable maximum weighted matching algorithm
    //! @sa mWeightedMatching
    inline void DisableWeightedMatching()
    {
        this->mWeightedMatching = 0;
    }

protected:
    //! @brief ... generate an error message from the error code
    //! @param error ... error code
    //! @return error message as std::string
    std::string GetErrorString(int error) const;

    //! @brief ... type of fill-in reducing odering of the coefficient matrix
    //! \sa setOrderingMETIS, setOrderingMinimumDegree
    /*!
        The attribute mOrderingType defines the type of fill-in reducing odering applied to the coefficient matrix.
        The minimum degree alogrithm (mOrderingType = 0) and the nested dissection alogrithm from the METIS-package
        (mOrderingType = 2) are available.
    */
    int mOrderingType;

    //! @brief ... maximum number of iterative refinement steps aplied in the solution phase
    //! @sa setNumRefinementSteps
    int mNumRefinementSteps;

    //! @brief ... absolute value of the order of magnitude of the perturbation for small pivots
    //! @sa setPivotingPerturbation
    /*!
        "This parameter instructs PARDISO how to handle small pivots or zero pivots for unsymmetric matrices
        and symmetric matrices." [Intel MKL manual]

        Small pivots will be perturbed with epsilon = 10^(-PivotingPerturbation).
        If the value of PivotingPerturbation is negative, default values will be used
        (perturbation of 1E-13 for unsymmetric matrices and 1E-8 for symmetric matrices).
    */
    int mPivotingPerturbation;

    //! @brief ... controls the scaling of the coefficient matrix
    //! @sa enableScaling, disableScaling
    /*!
         "PARDISO uses a maximum weight matching algorithm to permute large elements on the diagonal and to scale the
         matrix so that the diagonal elements are equal to 1 and the absolute values of the off-diagonal entries are
         less or equal to 1." [Intel MKL-Manual]

         If the value of scaling is negative, the scaling is only enabled for unsymmetric matrices (MKL default)
    */
    int mScaling;

    //! @brief ... controls maximum weighted matching algorithm
    //! @sa enableWeightedMatching, disableWeightedMatching
    /*!
        "PARDISO can use a maximum weighted matching algorithm to permute large elements close the diagonal. This strategy
        adds an additional level of reliability to our factorization ..." [Intel MKL-Manual]

        If the value of weightedMatching is negative, the maximum weighted matching algorithm is only enabled for
        unsymmetric matrices.
    */
    int mWeightedMatching;
#endif // HAVE_MKL_PARDISO
};
}
#endif // SPARSE_DIRECT_SOLVER_MKLPARDISO_H
