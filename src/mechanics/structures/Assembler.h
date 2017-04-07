#pragma once
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/constraints/Constraints.h"

namespace NuTo
{

class Assembler
{
public:
    Assembler();

#ifndef SWIG
    void BuildGlobalDofs(const std::vector<NodeBase*>& rNodes);

    //! @brief returns the number of constraint equations for a specific dof type
    //! @return number of constraints
    //! @param rDofType  dof type
    int ConstraintGetNumLinearConstraints(Node::eDof rDof) const;

#endif

    //! @brief returns the number of constraint equations for a specific dof type
    //! @return number of constraints
    //! @param rDofType  dof type
    int ConstraintGetNumLinearConstraints(std::string rDof) const;


    //! @brief returns the constraint vector after gauss elimination
    //! rConstraintMatrix*DOFS = RHS
    BlockFullVector<double> ConstraintGetRhsBeforeGaussElimination(double time) const;


    //! @brief calculates the right hand side of the constraint equations based on the mapping matrix and the rhs before the gauss elimination
    //! the result is stored internally in mConstraintRHS
    const BlockFullVector<double>& ConstraintGetRhsAfterGaussElimination() const
    {
        return mConstraintRhs;
    }

    void ConstraintUpdateRhs(double time);

    void AddEquation(NuTo::Node::eDof dof, Constraint::Equation equation);

    void AddEquations(NuTo::Node::eDof dof, std::vector<Constraint::Equation> equations);

    //! @brief summarizes information to dof numbering, active dof types, symmetric dof types, constant dof types
    DofStatus mDofStatus;
    
    //!brief ... renumbering of nodal DOFs required or not
    bool mNodeNumberingRequired;
    
    //! @brief constraint matrix relating the prescibed nodal unknowns to the free parameters
    BlockSparseMatrix mConstraintMatrix;

    //! @brief mapping matrix of the rhs to relate the rhs before the gauss elimination to the constraint matrix after
    // (mConstraintRHS (after elimination) = mConstraintMappingRHS *  mConstraintRHS (before elimination)
    // (the values of the RHS before elimination are stored at the individual constraints
    // the initial system is e.g.
    //[1 1 0]* [d1 d2 d3]^T = [rhs1]
    //[0 0 2]                 [rhs2]
    // this is replaced by
    //[1 1 0]* [d1 d2 d3]^T = rhs1 *[1] + rhs2 *[0]
    //[0 0 2]                       [0]         [1]
    // after gauss elimination and reordering this gives
    //[1 0 1]* [d1 d3 d2]^T = rhs1 *[1] + rhs2 *[0]
    //[0 1 0]                       [0]         [0.5]
    // as a consequence, the gauss elimination has only to be performed for a change of the constraint matrix
    // for a change of the rhs it is sufficient to recalculate the rhs from the above vectors
    // the mapping matrix [1,0; 0,0.5] is stored and the rhs is calculated from
    // mConstraintMappingRHS*mConstraintRHSBeforGaussElimination
    BlockSparseMatrix mConstraintMappingRhs;

    //! @brief right hand side of the constraint equations
    BlockFullVector<double> mConstraintRhs;

    Constraint::Constraints& GetConstraints()
    {
        return mConstraints;
    }

    const Constraint::Constraints& GetConstraints() const
    {
        return mConstraints;
    }

private:

    Constraint::Constraints mConstraints; 
};

} /* NuTo */
