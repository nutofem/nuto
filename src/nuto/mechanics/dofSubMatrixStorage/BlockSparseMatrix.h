#pragma once

#include <memory>

#include "nuto/mechanics/dofSubMatrixStorage/BlockStorageBase.h"
#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/math/SparseMatrixCSRGeneral_Def.h"
#include "nuto/math/SparseMatrixCSRSymmetric_Def.h"

namespace NuTo
{

//! @author Thomas Titscher, BAM
//! @date January 2016
//! @brief ... class for all block sparse matrices without any operations, except *BlockFullVector
class BlockSparseMatrix: public BlockStorageBase
{
public:
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    BlockSparseMatrix() {}
//    template<class Archive> void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief ctor
    //! @param rDofStatus ... reference to DofStatus for automatic matrix resizing
    //! @param rCanBeSymmetric ... true for e.g. M_JJ and M_KK, false for M_JK and M_KJ
    BlockSparseMatrix(const DofStatus& rDofStatus, bool rCanBeSymmetric = true);

    //! @brief copy ctor
    BlockSparseMatrix(const BlockSparseMatrix&  rOther);

#ifndef SWIG
    //! @brief move ctor
    //! @param rOther ... other BlockSparseMatrix
    BlockSparseMatrix(      BlockSparseMatrix&& rOther) = default;
#endif

    //! @brief allocates the submatrices based on the current dof configuration of the structure
    void AllocateSubmatrices();

    //! @brief set the right dimensions of the off-diagonal entries (and zero values)
    //! @remark In some cases, only the diagonal terms are set (constraint matrix). The off diagonals have to be resized for the operators to work properly.
    void FixOffDiagonalDimensions();

#ifndef SWIG
    //! @brief copy assignment
    //! @remark copies only active dof types
    BlockSparseMatrix& operator=(const BlockSparseMatrix&  rOther);

    //! @brief move assignment
    BlockSparseMatrix& operator=(      BlockSparseMatrix&& rOther) = default;

    friend std::ostream& operator<< (std::ostream &rOut, const NuTo::BlockSparseMatrix& rBlockSparseMatrix);

    //! @brief non-const access to the pair(rDofRow, rDofCol)
    //! @param rDofRow ... desired row
    //! @param rDofCol ... desired column
          SparseMatrixCSRVector2<double>& operator()(Node::eDof rDofRow, Node::eDof rDofCol);

    //! @brief const access to the pair(rDofRow, rDofCol)
    //! @param rDofRow ... desired row
    //! @param rDofCol ... desired column
    const SparseMatrixCSRVector2<double>& operator()(Node::eDof rDofRow, Node::eDof rDofCol) const;
#endif

    //! @brief adds the scaled rhs to this
    //! @param rRhs ... matrix to add
    //! @param rScalar ... scalar
    void AddScal(const BlockSparseMatrix& rRhs, double rScalar);

    //! @brief adds the scaled rhs as a diagonal matrix to this
    //! @param rRhs ... diagonal matrix (diagonal) to add
    //! @param rScalar ... scalar
    void AddScalDiag(const BlockFullVector<double>& rRhs, double rScalar);


    //! @brief allocates and returns a new BlockFullVector with the result
    //! @remark only calculates active dof values
    BlockFullVector<double> operator*(const BlockFullVector<double>& rRhs) const;

    //! @brief checks that each column of a row has the same number of rows, same vice versa
    void CheckDimensions() const;

#ifndef SWIG
    //! @brief gets the number of columns of the block storage for certain dof types
    //! @param rDofTypes ... dof types to count
    //! @return number of columns
    int GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const override;

    //! @brief gets the number of rows of the block storage for certain dof types
    //! @param rDofTypes ... dof types to count
    //! @return number of rows
    int GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const override;
#endif

    //! @brief prints submatrices and their dimensions
    void Info() const override;

    //! @brief returns true if the active dofs of the block matrix are symmetric
    //! @remark 1) only one dof type is active (and symmetric), 2) the diagonals are symmetric and the off-diagonals 0.
    bool HasSymmetricActiveDofs() const;

    //! @brief sets all sub matrices zero
    void SetZero();

    //! @brief inverts the matrix coefficient-wise
    void CwiseInvert();

#ifndef SWIG
    //! @brief resizes every member matrix according to the NumRowDofs and NumColumnDofs
    //! @param rNumRowDofsMap ... map containing the number of rows for each dof
    //! @param rNumColumnDofsMap ... map containing the number of columns for each dof
    void Resize(const std::map<Node::eDof, int>& rNumRowDofsMap, const std::map<Node::eDof, int>& rNumColumnDofsMap);
#endif

    //! @brief adds the GetNumEntries() for all active dof types
    int GetNumActiveEntires() const;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>  ExportToFullMatrix() const;
    NuTo::SparseMatrixCSRVector2General<double>               ExportToCSRVector2General() const;
    NuTo::SparseMatrixCSRGeneral<double>                      ExportToCSRGeneral() const;
    NuTo::SparseMatrixCSRVector2Symmetric<double>             ExportToCSRSymmetric() const;

    NuTo::SparseMatrixCSRVector2General<double> Get(std::string rDofRow, std::string rDofCol) const
    {
        auto& ref = (*this)(Node::DofToEnum(rDofRow), Node::DofToEnum(rDofCol));
        if (ref.IsSymmetric())
            return ref.AsSparseMatrixCSRVector2Symmetric(); // calls appropriate Vector2General ctor
        else
            return ref.AsSparseMatrixCSRVector2General();
    }

private:

    //! @brief storage using a unique_ptr
    //! @todo use std::make_unique as soon as the compiler version allows it
    std::map<std::pair<Node::eDof, Node::eDof>, std::unique_ptr<SparseMatrixCSRVector2<double>> > mData;

    bool mCanBeSymmetric;
};

} /* namespace NuTo */


