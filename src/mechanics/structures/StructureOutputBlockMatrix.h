#pragma once


#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/structures/StructureOutputBase.h"
#include <eigen3/Eigen/Sparse>


namespace NuTo
{
// Foward declarations
template <typename T>
class BlockFullVector;
template <typename T>
class BlockFullMatrix;
class ElementBase;

//! @author Volker Hirthammer
//! @date January 08, 2015
//! @brief ...
class StructureOutputBlockMatrix : public StructureOutputBase
{
private:
    using SparseMatrix = Eigen::SparseMatrix<double>;


    // Constructor / Destructor
    // ------------------------

public:
    //! @brief constructor
    //! @param rDofStatus: Reference to DofStatus needed for block matrices
    //! @param rAutomaticResize ... resizes the submatrices, if true
    StructureOutputBlockMatrix(const DofStatus& rDofStatus, bool rAutomaticResize = false);


    //! @brief copy constructor
    //! @remark StructureOutputBlockMatrix holds heavy data, no copies allowed
    StructureOutputBlockMatrix(const StructureOutputBlockMatrix& rOther) = default;


#ifndef SWIG

    //! @brief move constructor
    //! @param rOther ... other StructureOutputBlockMatrix
    StructureOutputBlockMatrix(StructureOutputBlockMatrix&& rOther) = default;

    //! @brief destructor
    ~StructureOutputBlockMatrix();


    // Operator overloads
    // ------------------


    //! @brief copy assignment operator
    StructureOutputBlockMatrix& operator=(const StructureOutputBlockMatrix& rOther) = default;


    //! @brief move assignment operator
    //! @param rOther ... other StructureOutputBlockMatrix
    StructureOutputBlockMatrix& operator=(StructureOutputBlockMatrix&& rOther) = default;

    StructureOutputBlockVector operator*(const StructureOutputBlockVector& rRhs) const;

    friend std::ostream& operator<<(std::ostream& rOut, const StructureOutputBlockMatrix& rStructureOutputBlockMatrix);


    // Member functions
    // ----------------

    void AddScal(const StructureOutputBlockMatrix& rOther, double rFactor);

    void AddScalDiag(const StructureOutputBlockVector& rOther, double rFactor);


    void AddElementMatrix(const ElementBase* rElementPtr, const NuTo::BlockFullMatrix<double>& rElementMatrix,
                          const NuTo::BlockFullVector<int>& rGlobalRowDofNumbers,
                          const NuTo::BlockFullVector<int>& rGlobalColumnDofNumbers, double mAddValueTolerance);

    void AddElementVectorDiagonal(const NuTo::BlockFullVector<double>& rElementVector,
                                  const NuTo::BlockFullVector<int>& rGlobalRowDofNumber, double mAddValueTolerance);


    //! @brief adds \f$(\boldsymbol{M}_{JJ} - \boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KJ} -
    //! \boldsymbol{M}_{JK}\,\boldsymbol{C}_{mat} +
    //! \boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KK}\,\boldsymbol{C}_{mat})\,c\f$ to rHessian
    //! @remark only calculates active dof types
    //! @param rHessian ... global hessian
    //! @param rCmat ... constraint matrix
    //! @param rScalar ... option to scale the terms
    void ApplyCMatrixScal(BlockSparseMatrix& rHessian, const BlockSparseMatrix& rCmat, double rScalar) const;

    void ApplyCMatrix(const BlockSparseMatrix& rCmat);


    //! @brief resizes every member matrix according to numActiveDofs/numDependentDofs
    void Resize(const std::map<Node::eDof, int>& rNumActiveDofsMap,
                const std::map<Node::eDof, int>& rNumDependentDofsMap);

    //! @brief checks the dimensions of all the submatrices for consistency
    //! @return true if matrix dimensions are consistent, false if not
    void CheckDimensions() const;

    StructureOutputBlockMatrix& AsStructureOutputBlockMatrix() override
    {
        return *this;
    }

    //! @brief inverts the matrix coefficient-wise
    void CwiseInvert();

    void SetZero() override;

    SparseMatrix ExportToEigenSparseMatrix() const;


#endif

private:
    void InsertSubMatrix(SparseMatrix& rMat, const SparseMatrix& subMat, const int startRowId,
                         const int startColId) const;
    // Member variables
    // ----------------

public:
    NuTo::BlockSparseMatrix JJ;
    NuTo::BlockSparseMatrix JK;
    NuTo::BlockSparseMatrix KJ;
    NuTo::BlockSparseMatrix KK;
};


} // namespace NuTo
