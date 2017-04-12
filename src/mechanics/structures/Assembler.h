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

    //! @brief ctor
    Assembler();

    //! @brief builds the global dof numbering depending on the constraints
    //! sets the members mConstraintMatrix, mConstraintMappingRhs and mConstraintRhs [for time t=0]
    //! @param nodes all the nodes included in the global dof numbering
    void BuildGlobalDofs(const std::vector<NodeBase*>& nodes);

    //! @brief getter for mConstraintRhs, set via ConstraintUpdateRhs
    const BlockFullVector<double>& GetConstraintRhs() const
    {
        return mConstraintRhs;
    }

    //! @brief getter for mConstraintMatrix, set via BuildGlobalDofs
    const BlockSparseMatrix& GetConstraintMatrix() const
    {
        return mConstraintMatrix;
    }
    
    //! @brief calculates the right hand side of the constraint equations based on the mapping matrix and the rhs before the gauss elimination
    //! the result is stored internally in mConstraintRHS
    //! @param time global time
    void ConstraintUpdateRhs(double time);

    //! @brief throws if the nodes or constraints equations have changed
    void ThrowIfRenumberingRequred() const;

    //! @brief returns true if a node renumbering is required (nodes or constraints changed)
    bool RenumberingRequired() const
    {
        return (mConstraints.HaveChanged() or mNodeVectorChanged);
    }

    //! @brief sets the state variable mNodeVectorChanged to true to indicate
    //! that a node renumbering is required
    void SetNodeVectorChanged()
    {
        mNodeVectorChanged = true;
    }

    //! @brief non-const getter for mConstraints
    Constraint::Constraints& GetConstraints()
    {
        return mConstraints;
    }

    //! @brief getter for mConstraints
    const Constraint::Constraints& GetConstraints() const
    {
        return mConstraints;
    }

    //! @brief summarizes information to dof numbering, active dof types, symmetric dof types, constant dof types
    //! @TODO: should not be pulbic
    DofStatus mDofStatus;

private:
    
    //! @brief builds the constraint rhs vector before the gauss elimination evaluated at time
    //! @param time global time
    //! @return constraint rhs before gauss elimination
    BlockFullVector<double> BuildRhsBeforeGaussElimination(double time) const;

    //! @brief stores the constraints
    Constraint::Constraints mConstraints; 
    
    //!brief ... renumbering of nodal DOFs required or not
    bool mNodeVectorChanged = false;
    
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
};

} /* NuTo */
