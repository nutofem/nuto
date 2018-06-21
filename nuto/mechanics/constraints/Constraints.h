#pragma once
#include <set>
#include <eigen3/Eigen/Sparse>
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/constraints/Equation.h"

namespace NuTo
{
namespace Constraint
{

//! @brief represents a numbering where first the independent dofs are consecutively numbered
//! leaving out the dependent ones and then the dependent dof numbers are appended.
//!
//! mNumJ number of independent Dofs
//! mNumK number of dependent Dofs
//!
//! The order of the dependent dofs is given by the order of the constraint equations
//!
//! Example:
//!
//! dofNumber                   0  1  2*  3  4  5*  6  (* = constrained)
//! constraintEqNr                    1         0
//!
//! JKNumbering
//! (independent dofs)          0  1      2  3      4
//!
//! (dependent dofs)
//! appended and ordered
//! like constraint equations         6         5
//!
struct JKNumbering
{
    JKNumbering(Eigen::VectorXi indices, int numK)
        : mIndices(indices)
        , mNumJ(indices.size() - numK)
        , mNumK(numK)
    {
    }

    Eigen::VectorXi mIndices;
    int mNumJ;
    int mNumK;
};

//! @brief stores constraint equations, separated by their dof type
class Constraints
{
    using Equations = std::vector<Equation>;

public:
    //! @brief adds an equation
    //! @param dof dof type
    //! @param equation constraint equation
    void Add(DofType dof, Equation equation);

    //! @brief adds multiple equations
    //! @param dof dof type
    //! @param equations constraint equations
    void Add(DofType dof, std::vector<Equation> equations);

    //! @brief builds the time dependent constraint rhs vector for a specific dof type
    //! @param dof dof type
    //! @param time global time
    //! @return constraint rhs vector
    Eigen::VectorXd GetRhs(DofType dof, double time) const;

    //! @brief builds the rhs of the constraint equations in a sparse global vector
    //! @param dof dof type
    //! @param numDofs total number of dofs for a specfic dof type
    //! @param time time at which the rhs vector should be evaluated
    //! @return sparse rhs vector of the constraint equations rhs(time)
    Eigen::SparseVector<double> GetSparseGlobalRhs(DofType dof, int numDofs, double time) const;

    //! @brief builds a numbering where first the independent dofs are consecutively numbered
    //! leaving out the dependent ones and then the dependent dof numbers are appended. The
    //! order of the dependent dofs is given by the order of the constraint equations
    //!
    //! This numbering is used in BuildUnitConstraintMatrix
    //! @param dof dof type
    //! @param numDofs number dofs for the dof type
    //! @return index vector with first independent dof numbers and then dependent
    JKNumbering GetJKNumbering(DofType dof, int numDofs) const;

    //! @brief builds a sparse matrix containing the constraint terms for a specific dof type and a unit matrix for the
    //! independent dofs
    //! @param dof dof type
    //! @param numDofs number dofs for the dof type, required for a proper resize of the sparse matrix
    //! @return sparse matrix containing the constraint terms
    Eigen::SparseMatrix<double> BuildUnitConstraintMatrix(DofType dof, int numDofs) const;

    //! @brief builds a sparse matrix to extract the independent values from the full vector of dofs
    //! @param dof dof type
    //! @param numDofs number dofs for the dof type, required for a proper resize of the sparse matrix
    //! @return sparse matrix
    Eigen::SparseMatrix<double> BuildUnitConstraintMatrixInv(DofType dof, int numDofs) const;

    //! @brief calculates the number of constraint equations for a specific dof type
    //! @param dof dof type
    //! @return number of constraint equations
    int GetNumEquations(DofType dof) const;

    //!@brief gets the specified equation
    //! @param dof doftype of the equation
    //! @param equationNumber number of the equation
    //! @return selected equation
    const Equation& GetEquation(DofType dof, int equationNumber) const;

private:
    //! @brief dof-wise storage of constraint equations
    DofContainer<Equations> mEquations;


    //! @brief stores the terms of all existing equations in ordered containers
    //! @remark A naive check of all existing terms in mEquations is O(N^2) since we have to traverse them whenever an
    //! equation is added. Adding >10000 equations this way is not feasible. This class performs those checks in O(N log
    //! N) and everything is fine again.
    class TermChecker
    {
    public:
        //! @brief Check if the `e` collides with existing terms.
        //! @remark throws if
        //          - the dependent term of `e` is in any existing equation
        //          - any term of `e` is in the depenend terms of any existing equation
        void CheckEquation(Equation e);

    private:
        struct TermCompare
        {
            bool operator()(const Term& lhs, const Term& rhs) const;
        };

        std::set<Term, TermCompare> mDependentTerms;
        std::set<Term, TermCompare> mIndependentTerms;
    } mTermChecker;
};


} /* Constaint */
} /* NuTo */
