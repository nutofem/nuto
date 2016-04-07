#pragma once

#include "nuto/mechanics/dofSubMatrixStorage/BlockStorageBase.h"
#include "nuto/math/FullMatrix_Def.h"

namespace NuTo
{
template <typename T> class BlockFullVector;
//! @author Thomas Titscher, BAM
//! @date January 2016
//! @brief ... class for all block full matrices, only for storing data, no calculations
template <typename T>
class BlockFullMatrix: public BlockStorageBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    BlockFullMatrix() {}
//    template<class Archive> void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

public:

    //! @brief ctor
    //! @param rDofStatus ... reference to DofStatus for automatic matrix resizing
    BlockFullMatrix(const DofStatus& rDofStatus);

    //! @brief copy ctor
    BlockFullMatrix(const BlockFullMatrix& rOther) = default;

#ifndef SWIG
    //! @brief move ctor
    BlockFullMatrix(BlockFullMatrix&& rOther) = default;

    //! @brief non-const access
          NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& operator()(Node::eDof rDofRow, Node::eDof rDofCol);

    //! @brief const access
    const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& operator()(Node::eDof rDofRow, Node::eDof rDofCol) const;

    //! @brief move assignment
    BlockFullMatrix& operator =(BlockFullMatrix&& rOther) = default;


    //! @brief move assignment
    BlockFullMatrix& operator =(const BlockFullMatrix& rOther) = default;

    //! @brief operator +=
    //! @remark only modifies active dof types
    BlockFullMatrix& operator+=(const BlockFullMatrix& rRhs);

    //! @brief operator -=
    //! @remark only modifies active dof types
    BlockFullMatrix& operator-=(const BlockFullMatrix& rRhs);

    friend NuTo::BlockFullMatrix<T> operator+(NuTo::BlockFullMatrix<T> rLhs, const NuTo::BlockFullMatrix<T>& rRhs) { return std::move(rLhs += rRhs); }
    friend NuTo::BlockFullMatrix<T> operator-(NuTo::BlockFullMatrix<T> rLhs, const NuTo::BlockFullMatrix<T>& rRhs) { return std::move(rLhs -= rRhs); }

    template <typename T2>
    friend std::ostream& operator<< (std::ostream &rOut, const NuTo::BlockFullMatrix<T2>& rBlockVector);
#endif

    //! @brief prints submatrices and their dimensions
    void Info() const override;

    //! @brief checks that each column of a row has the same number of rows, same vice versa
    void CheckDimensions() const;

#ifndef SWIG

    //! @brief gets the number of columns of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of columns
    int GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const override;

    //! @brief gets the number of rows of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of rows
    int GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const override;


#endif // SWIG

    NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> Export() const;

    NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> Get(std::string rDofRow, std::string rDofCol) const
    {
        return (*this)(Node::DofToEnum(rDofRow), Node::DofToEnum(rDofCol));
    }

private:
    std::map<std::pair<Node::eDof, Node::eDof>, NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>> mData;
};


} /* namespace NuTo */
