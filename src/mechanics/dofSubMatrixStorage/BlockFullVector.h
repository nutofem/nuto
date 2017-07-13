#pragma once


#include "mechanics/dofSubMatrixStorage/BlockStorageBase.h"
#include "mechanics/nodes/DofHash.h"
#include "eigen3/Eigen/Dense"
#include <ostream>
#include <unordered_map>
#include <map>

namespace NuTo
{
class BlockScalar;
class SerializeStreamOut;
class SerializeStreamIn;

//! @author Thomas Titscher, BAM
//! @date January 2016
//! @brief ... class for all block vectors with basic operators +,-,*(scalar)
template <typename T>
class BlockFullVector : public BlockStorageBase
{

public:
    //! @brief ctor
    //! @param rDofStatus ... reference to DofStatus for automatic matrix resizing
    BlockFullVector(const DofStatus& rDofStatus);

    //! @brief copy constructor
    BlockFullVector(const BlockFullVector& rOther);

    //! @brief destructor
    ~BlockFullVector();

#ifndef SWIG
    //! @brief move constructor
    BlockFullVector(BlockFullVector&& rOther);
#endif

    //! @brief import constructor
    //! @param rData ... vector to import from
    //! @param rDofStatus ... reference to DofStatus for automatic matrix resizing
    //! @param rAreActiveDofValues ... true if the rData represents a vector of active dof values
    BlockFullVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& rData, const DofStatus& rDofStatus,
                    bool rAreActiveDofValues = true);

    //! @brief allocates the subvectors based on the current dof configuration of the structure
    void AllocateSubvectors();


#ifndef SWIG
    //! @brief copy assignment
    //! @remark only copies the active dof types
    BlockFullVector& operator=(const BlockFullVector& rOther);

    //! @brief move assignment
    //! @remark moves all values, moving only the active dof types is somehow slower.
    BlockFullVector& operator=(BlockFullVector&& rOther);

    //! @brief non-const access
    Eigen::Matrix<T, Eigen::Dynamic, 1>& operator[](Node::eDof rDofRow);

    //! @brief const access
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& operator[](Node::eDof rDofRow) const;


    //! @brief operator +=
    //! @remark only modifies active dof types
    BlockFullVector& operator+=(const BlockFullVector& rRhs);

    //! @brief operator -=
    //! @remark only modifies active dof types
    BlockFullVector& operator-=(const BlockFullVector& rRhs);

    //! @brief operator *=
    //! @remark only modifies active dof types
    BlockFullVector& operator*=(double rScalar);

    //! @brief operator *=
    //! @remark only modifies active dof types
    BlockFullVector& operator/=(double rScalar);


    friend NuTo::BlockFullVector<T> operator+(NuTo::BlockFullVector<T> rLhs, const NuTo::BlockFullVector<T>& rRhs)
    {
        return std::move(rLhs += rRhs);
    }
    friend NuTo::BlockFullVector<T> operator-(NuTo::BlockFullVector<T> rLhs, const NuTo::BlockFullVector<T>& rRhs)
    {
        return std::move(rLhs -= rRhs);
    }
    friend NuTo::BlockFullVector<T> operator*(NuTo::BlockFullVector<T> rLhs, double rScalar)
    {
        return std::move(rLhs *= rScalar);
    }
    friend NuTo::BlockFullVector<T> operator*(double rScalar, NuTo::BlockFullVector<T> rRhs)
    {
        return std::move(rRhs *= rScalar);
    }
    friend NuTo::BlockFullVector<T> operator/(NuTo::BlockFullVector<T> rLhs, double rScalar)
    {
        return std::move(rLhs /= rScalar);
    }
    friend NuTo::BlockFullVector<T> operator/(double rScalar, NuTo::BlockFullVector<T> rRhs)
    {
        return std::move(rRhs /= rScalar);
    }

    template <typename T2>
    friend std::ostream& operator<<(std::ostream& rOut, const NuTo::BlockFullVector<T2>& rBlockVector);

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    void NuToSerializeSave(SerializeStreamOut& rStream)
    {
        SerializeBlockFullVector(rStream);
    }

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    void NuToSerializeLoad(SerializeStreamIn& rStream)
    {
        SerializeBlockFullVector(rStream);
    }

#endif

    //! @brief comparision, checks equality of all sub vectors
    inline bool operator==(const BlockFullVector& rOther)
    {
        return mData == rOther.mData;
    }

    //! @brief comparision, checks !equality of all sub vectors
    inline bool operator!=(const BlockFullVector& rOther)
    {
        return !(*this == rOther);
    }

    //! @brief Imports the active dof type values from a vector
    void Import(const Eigen::Matrix<T, Eigen::Dynamic, 1>& rToImport);

    //! @brief prints subvectors and their dimensions
    void Info() const override;


#ifndef SWIG
    //! @brief gets the number of columns of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of columns
    int GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const override;

    //! @brief gets the number of rows of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of rows
    int GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const override;

    //! @brief resizes every member vector according to the NumRowDofs
    //! @param rNumRowDofsMap ... map containing the number of rows for each dof
    void Resize(const std::map<Node::eDof, int>& rNumRowDofsMap);

    //! @brief Calculates the L2 norm of the block vector for each dof
    BlockScalar CalculateNormL2();

    //! @brief Calculates the infinity norm of the block vector for each dof
    BlockScalar CalculateInfNorm() const;

#endif

    //! @brief sets every sub vector to zero
    void SetZero();


    //! @brief Exports the active dof type values to a vector
    Eigen::Matrix<T, Eigen::Dynamic, 1> Export() const;

    //! @brief Exports a specific dof type values to a vector
    //! @param rDofRow specific dof type string
    Eigen::Matrix<T, Eigen::Dynamic, 1> Get(std::string rDofRow) const;

private:
    //! @brief defines the serialization of this class
    //! @param rStream serialize input/output stream
    template <typename TStream>
    void SerializeBlockFullVector(TStream& rStream);

    std::unordered_map<Node::eDof, Eigen::Matrix<T, Eigen::Dynamic, 1>, Node::eDofHash> mData;
};

} /* namespace NuTo */
