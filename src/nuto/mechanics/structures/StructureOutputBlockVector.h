#pragma once

#include <iostream>
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/structures/StructureOutputBase.h"



namespace NuTo
{
    // Foward declarations
    class BlockSparseMatrix;




//! @author Volker Hirthammer
//! @date January 25, 2015
//! @brief ...
class StructureOutputBlockVector: public StructureOutputBase
{

public:

    //! @brief constructor
    //! @param rDofStatus: Reference to DofStatus needed for block vectors
    //! @param rAutomaticResize ... resizes the submatrices, if true
    StructureOutputBlockVector(const DofStatus& rDofStatus, bool rAutomaticResize = false);

    //! @brief copy constructor
    StructureOutputBlockVector(const StructureOutputBlockVector&  rOther) = default;

#ifndef SWIG
    //! @brief move constructor
    //! @param rOther ... other StructureOutputBlockVector
    StructureOutputBlockVector(      StructureOutputBlockVector&& rOther) = default;

    //! @brief destructor
    ~StructureOutputBlockVector() = default;



    // Operator overloads
    // ------------------


    //! @brief copy assignment operator
    StructureOutputBlockVector& operator=(const StructureOutputBlockVector&  rOther) = default;


    //! @brief move assignment operator
    //! @param rOther ... other StructureOutputBlockVector
    StructureOutputBlockVector& operator=(      StructureOutputBlockVector&& rOther) = default;

    StructureOutputBlockVector& operator+= (const StructureOutputBlockVector&  rRhs);
    StructureOutputBlockVector& operator-= (const StructureOutputBlockVector&  rRhs);
    StructureOutputBlockVector& operator*= (double rRhs);
    StructureOutputBlockVector& operator/= (double rRhs);


    friend StructureOutputBlockVector operator+ (StructureOutputBlockVector rLhs, const StructureOutputBlockVector& rRhs)  { return std::move(rLhs += rRhs); }
    friend StructureOutputBlockVector operator- (StructureOutputBlockVector rLhs, const StructureOutputBlockVector& rRhs)  { return std::move(rLhs -= rRhs); }
    friend StructureOutputBlockVector operator* (StructureOutputBlockVector rLhs, double rRhs)
    {
        return std::move(rLhs *= rRhs);
    }
    friend StructureOutputBlockVector operator* (double rLhs, StructureOutputBlockVector rRhs)
    {
        return std::move(rRhs *= rLhs);
    }
    friend StructureOutputBlockVector operator/ (StructureOutputBlockVector rLhs, double rRhs)
    {
        return std::move(rLhs /= rRhs);
    }
    friend StructureOutputBlockVector operator/ (double rLhs, StructureOutputBlockVector rRhs)
    {
        return std::move(rRhs /= rLhs);
    }

    friend std::ostream& operator<< (std::ostream &rOut, const StructureOutputBlockVector& rStructureOutputBlockVector);

    // Member functions
    // ----------------

    void AddElementVector(
            const NuTo::BlockFullVector<double>& rElementVector,
            const NuTo::BlockFullVector<int>& rGlobalRowDofNumbers);

    //! @brief adds \f$(\boldsymbol{R}_{J} - \boldsymbol{C}_{mat}^T\,\boldsymbol{R}_{K}),c\f$ to rResidual
    //! @remark only calculates active dof types
    //! @param rResidual ... global residual
    //! @param rCmat ... constraint matrix
    void ApplyCMatrix(BlockFullVector<double>& rResidual, const BlockSparseMatrix& rCmat) const;

    //! @brief adds \f$(- \boldsymbol{C}_{mat}^T\,\boldsymbol{R}_{K}),c\f$ to this.J
    //! @remark only calculates active dof types
    //! @param rCmat ... constraint matrix
    void ApplyCMatrix(const BlockSparseMatrix& rCmat);

    //! @brief resizes every member vector according to numActiveDofs/numDependentDofs
    void Resize(const std::map<Node::eDof, int>& rNumActiveDofsMap, const std::map<Node::eDof, int>& rNumDependentDofsMap);

    //! @brief sets each member vector to zero
    void SetZero() override;

    StructureOutputBlockVector& AsStructureOutputBlockVector() override
    {
        return *this;
    }

#endif // SWIG

public:

    NuTo::BlockFullVector<double> J;
    NuTo::BlockFullVector<double> K;
};


} // namespace NuTo
