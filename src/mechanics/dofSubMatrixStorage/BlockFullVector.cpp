#include <iostream>

#include "base/Exception.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/nodes/NodeEnum.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"

template <typename T>
NuTo::BlockFullVector<T>::BlockFullVector(const DofStatus& rDofStatus)
    : BlockStorageBase(rDofStatus)
{
    AllocateSubvectors();
}

template <typename T>
NuTo::BlockFullVector<T>::~BlockFullVector()
{
}

template <typename T>
NuTo::BlockFullVector<T>::BlockFullVector(const NuTo::BlockFullVector<T>& rOther)
    : BlockStorageBase(rOther.mDofStatus)
    , mData(rOther.mData)
{
}

template <typename T>
NuTo::BlockFullVector<T>::BlockFullVector(NuTo::BlockFullVector<T>&& rOther)
    : BlockStorageBase(rOther.mDofStatus)
    , mData(std::move(rOther.mData))
{
}

template <typename T>
NuTo::BlockFullVector<T>::BlockFullVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& rData, const DofStatus& rDofStatus,
                                          bool rAreActiveDofValues)
    : BlockStorageBase(rDofStatus)
{
    AllocateSubvectors();
    if (rAreActiveDofValues)
        Resize(rDofStatus.GetNumActiveDofsMap());
    else
        Resize(rDofStatus.GetNumDependentDofsMap());

    Import(rData);
}

template <typename T>
void NuTo::BlockFullVector<T>::AllocateSubvectors()
{
    mData.clear();
    for (auto dof : mDofStatus.GetDofTypes())
    {
        mData[dof] = Eigen::Matrix<T, Eigen::Dynamic, 1>();
    }
}

template <typename T>
NuTo::BlockFullVector<T>& NuTo::BlockFullVector<T>::operator=(const BlockFullVector<T>& rOther)
{
    for (auto dof : mDofStatus.GetActiveDofTypes())
        mData[dof] = rOther[dof];
    return *this;
}

template <typename T>
NuTo::BlockFullVector<T>& NuTo::BlockFullVector<T>::operator=(NuTo::BlockFullVector<T>&& rOther)
{
    mData = std::move(rOther.mData);
    return *this;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>& NuTo::BlockFullVector<T>::operator[](Node::eDof rDofRow)
{
    auto data = mData.find(rDofRow);
    assert(data != mData.end());
    return (*data).second;
}


template <typename T>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& NuTo::BlockFullVector<T>::operator[](Node::eDof rDofRow) const
{
    auto data = mData.find(rDofRow);
    assert(data != mData.end());
    return (*data).second;
}


template <typename T>
NuTo::BlockFullVector<T>& NuTo::BlockFullVector<T>::operator+=(const BlockFullVector<T>& rRhs)
{
    for (auto dof : mDofStatus.GetActiveDofTypes())
        (*this)[dof] += rRhs[dof];
    return *this;
}

template <typename T>
NuTo::BlockFullVector<T>& NuTo::BlockFullVector<T>::operator-=(const BlockFullVector<T>& rRhs)
{
    for (auto dof : mDofStatus.GetActiveDofTypes())
        (*this)[dof] -= rRhs[dof];
    return *this;
}

template <typename T>
NuTo::BlockFullVector<T>& NuTo::BlockFullVector<T>::operator*=(double rScalar)
{
    for (auto dof : mDofStatus.GetActiveDofTypes())
        (*this)[dof] *= rScalar;
    return *this;
}

template <typename T>
NuTo::BlockFullVector<T>& NuTo::BlockFullVector<T>::operator/=(double rScalar)
{
    for (auto dof : mDofStatus.GetActiveDofTypes())
        (*this)[dof] /= rScalar;
    return *this;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> NuTo::BlockFullVector<T>::Export() const
{
    Eigen::Matrix<T, Eigen::Dynamic, 1> result(GetNumActiveRows());

    int blockStartIndex = 0;
    for (auto dof : mDofStatus.GetActiveDofTypes())
    {
        int numRows = (*this)[dof].rows();
        result.block(blockStartIndex, 0, numRows, 1) = (*this)[dof];
        blockStartIndex += numRows;
    }
    return result;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> NuTo::BlockFullVector<T>::Get(std::string rDofRow) const
{
    return (*this)[Node::DofToEnum(rDofRow)];
}

template <typename T>
void NuTo::BlockFullVector<T>::Import(const Eigen::Matrix<T, Eigen::Dynamic, 1>& rToImport)
{
    if (GetNumActiveRows() != rToImport.rows())
    {
        this->Info();
        throw NuTo::Exception(std::string("[") + __PRETTY_FUNCTION__ +
                                       "] BlockFullVector must be sized to the right dimensions");
    }

    int blockStartIndex = 0;
    for (auto dof : mDofStatus.GetActiveDofTypes())
    {
        int numRows = (*this)[dof].rows();
        (*this)[dof] = rToImport.block(blockStartIndex, 0, numRows, 1);
        blockStartIndex += numRows;
    }
}


template <typename T>
void NuTo::BlockFullVector<T>::Info() const
{
    int minLength = 30;
    std::cout << "Num sub vectors: " << mData.size() << std::endl;
    for (auto& pair : mData)
    {
        std::string dofTypes = "[ " + Node::DofToString(pair.first) + " ]:";
        int numAdditionalBlanks = std::max(0, minLength - (int)dofTypes.length());
        const std::string& additionalBlanks = std::string(numAdditionalBlanks, ' ');
        const auto& vector = pair.second;

        std::cout << dofTypes << additionalBlanks << "(" << vector.rows() << "x" << vector.cols() << ")" << std::endl;
    }
}


template <typename T>
int NuTo::BlockFullVector<T>::GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const
{
    return 1;
}


template <typename T>
int NuTo::BlockFullVector<T>::GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const
{
    int numRows = 0;
    for (auto dof : rDofTypes)
    {
        numRows += (*this)[dof].rows();
    }
    return numRows;
}

template <typename T>
void NuTo::BlockFullVector<T>::Resize(const std::map<Node::eDof, int>& rNumRowDofsMap)
{
    const auto& dofTypes = mDofStatus.GetDofTypes();
    for (auto dofRow : dofTypes)
    {
        (*this)[dofRow].resize(rNumRowDofsMap.at(dofRow));
        (*this)[dofRow].setZero();
    }
}

template <typename T>
void NuTo::BlockFullVector<T>::SetZero()
{
    for (auto& pair : mData)
        pair.second.setZero();
}

//! @brief Calculates the L2 norm of the block vector for each dof
template <typename T>
NuTo::BlockScalar NuTo::BlockFullVector<T>::CalculateNormL2()
{
    BlockScalar dofWiseNorm(mDofStatus);
    for (auto dof : mDofStatus.GetActiveDofTypes())
    {
        dofWiseNorm[dof] = mData[dof].norm();
    }
    return dofWiseNorm;
}

template <typename T>
NuTo::BlockScalar NuTo::BlockFullVector<T>::CalculateInfNorm() const
{
    BlockScalar dofWiseNorm(mDofStatus);
    for (auto dof : mDofStatus.GetActiveDofTypes())
    {
        dofWiseNorm[dof] = (*this)[dof].cwiseAbs().maxCoeff();
    }
    return dofWiseNorm;
}

template <typename T>
template <typename TStream>
void NuTo::BlockFullVector<T>::SerializeBlockFullVector(TStream& rStream)
{
    for (auto& data : mData)
    {
        rStream.Serialize(data.second);
    }
}

namespace NuTo
{
////! @brief stream operator for outputs with cout or files
template <typename T>
std::ostream& operator<<(std::ostream& rOut, const BlockFullVector<T>& rBlockVector)
{
    for (auto dof : rBlockVector.mDofStatus.GetActiveDofTypes())
    {
        rOut << "[" << Node::DofToString(dof) << "]" << std::endl;
        rOut << rBlockVector[dof] << std::endl;
    }
    return rOut;
}
} // namespace NuTo

template std::ostream& NuTo::operator<<(std::ostream& rOut, const NuTo::BlockFullVector<double>& rBlockVector);
template std::ostream& NuTo::operator<<(std::ostream& rOut, const NuTo::BlockFullVector<int>& rBlockVector);

template class NuTo::BlockFullVector<double>;
template class NuTo::BlockFullVector<int>;

template void NuTo::BlockFullVector<int>::SerializeBlockFullVector(NuTo::SerializeStreamIn&);
template void NuTo::BlockFullVector<int>::SerializeBlockFullVector(NuTo::SerializeStreamOut&);
template void NuTo::BlockFullVector<double>::SerializeBlockFullVector(NuTo::SerializeStreamIn&);
template void NuTo::BlockFullVector<double>::SerializeBlockFullVector(NuTo::SerializeStreamOut&);
