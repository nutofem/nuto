#pragma once
#include <boost/ptr_container/ptr_map.hpp>
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/constraints/ConstraintBase.h"

namespace NuTo
{

class Assembler
{
public:
    Assembler(const DofStatus& rDofStatus)
        : mConstraintMatrix(rDofStatus, false)
        , mConstraintMappingRHS(rDofStatus, false)
        , mConstraintRHS(rDofStatus)
    {
    }

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
    BlockSparseMatrix mConstraintMappingRHS;

    //! @brief right hand side of the constraint equations
    BlockFullVector<double> mConstraintRHS;


    //! @brief ... map storing the constraints
    //! @sa ConstraintBase
    boost::ptr_map<int,ConstraintBase> mConstraintMap;

};

} /* NuTo */
