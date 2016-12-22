#include <sstream>
#include "mechanics/MechanicsException.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/nodes/NodeEnum.h"
#include "math/FullMatrix.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>
#endif // ENABLE_SERIALIZATION


template<typename T>
NuTo::BlockFullMatrix<T>::BlockFullMatrix(const DofStatus& rDofStatus) : BlockStorageBase(rDofStatus)
{
    const auto& dofTypes = mDofStatus.GetDofTypes();
    for (auto dofRow : dofTypes)
        for (auto dofCol : dofTypes)
            mData[std::make_pair(dofRow, dofCol)]; // access once to allocate
}

template<typename T>
NuTo::BlockFullMatrix<T>::~BlockFullMatrix()
{}


template<typename T>
NuTo::BlockFullMatrix<T>::BlockFullMatrix(const NuTo::BlockFullMatrix<T> &rOther)
    : BlockStorageBase(rOther.mDofStatus),
      mData(rOther.mData)
{}

template<typename T>
NuTo::BlockFullMatrix<T>::BlockFullMatrix(NuTo::BlockFullMatrix<T>&& rOther)
    : BlockStorageBase(rOther.mDofStatus),
      mData(std::move(rOther.mData))
{}

template<typename T>
NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& NuTo::BlockFullMatrix<T>::operator ()(Node::eDof rDofRow, Node::eDof rDofCol)
{
    auto data = mData.find(std::make_pair(rDofRow, rDofCol));
    assert (data != mData.end());
    return (*data).second;
}

template<typename T>
const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& NuTo::BlockFullMatrix<T>::operator ()(Node::eDof rDofRow, Node::eDof rDofCol) const
{
    auto data = mData.find(std::make_pair(rDofRow, rDofCol));
    assert (data != mData.end());
    return (*data).second;
}

template<typename T>
NuTo::BlockFullMatrix<T>& NuTo::BlockFullMatrix<T>::operator =(const NuTo::BlockFullMatrix<T> &rOther)
{
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    for (auto dofRow : activeDofTypes)
        for (auto dofCol : activeDofTypes)
            (*this)(dofRow, dofCol) = rOther(dofRow, dofCol);
    return *this;
}

template<typename T>
NuTo::BlockFullMatrix<T>& NuTo::BlockFullMatrix<T>::operator =(NuTo::BlockFullMatrix<T>&& rOther)
{
    mData = std::move(rOther.mData);
    return *this;
}

template<typename T>
NuTo::BlockFullMatrix<T>& NuTo::BlockFullMatrix<T>::operator+=(const BlockFullMatrix& rRhs)
{
    for (auto dofRow : mDofStatus.GetActiveDofTypes())
        for (auto dofCol : mDofStatus.GetActiveDofTypes())
            (*this)(dofRow, dofCol) += rRhs(dofRow, dofCol);
    return *this;
}

template<typename T>
NuTo::BlockFullMatrix<T>& NuTo::BlockFullMatrix<T>::operator-=(const BlockFullMatrix& rRhs)
{
    for (auto dofRow : mDofStatus.GetActiveDofTypes())
        for (auto dofCol : mDofStatus.GetActiveDofTypes())
            (*this)(dofRow, dofCol) -= rRhs(dofRow, dofCol);
    return *this;
}



template<typename T>
inline void NuTo::BlockFullMatrix<T>::Info() const
{
    int minLength = 60;
    std::cout << "Num sub matrices: " << mData.size() << std::endl;
    for (auto& pair : mData)
    {
        std::string dofTypes = "[ " + Node::DofToString(pair.first.first) + " , " + Node::DofToString(pair.first.second) + " ]:";
        int numAdditionalBlanks = std::max(0, minLength - (int)dofTypes.length());
        const std::string& additionalBlanks = std::string(numAdditionalBlanks, ' ');
        const auto& matrix = pair.second;

        std::cout << dofTypes << additionalBlanks << "(" << matrix.GetNumRows() << "x" << matrix.GetNumColumns() << ")" << std::endl;
    }
}

template<typename T>
void NuTo::BlockFullMatrix<T>::CheckDimensions() const
{
    const auto& dofTypes = mDofStatus.GetDofTypes();
    /*
     * check for the same numRows in each row
     */
    for (auto dofRow : dofTypes)
    {
        auto numRowsReference = (*this)(dofRow, dofRow).GetNumRows();
        for (auto dofCol : dofTypes)
        {
            int numRows = (*this)(dofRow, dofCol).GetNumRows();
            if (numRows != numRowsReference)
            {
                std::stringstream s;
                s << "[" << __PRETTY_FUNCTION__ << "] Submatrix row dimension mismatch. \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofCol) << ") has " << numRows << " rows \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofRow) << ") has " << numRowsReference << " rows \n";
                throw MechanicsException(s.str());
            }
        }
    }

    /*
     * check for the same numCols in each column
     */
    for (auto dofCol : dofTypes)
    {
        auto numColsReference = (*this)(dofCol, dofCol).GetNumColumns();
        for (auto dofRow : dofTypes)
        {
            int numCols = (*this)(dofRow, dofCol).GetNumColumns();
            if (numCols != numColsReference)
            {
                std::stringstream s;
                s << "[" << __PRETTY_FUNCTION__ << "] Submatrix column dimension mismatch. \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofCol) << ") has " << numCols << " columns \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofRow) << ") has " << numColsReference << " columns \n";
                throw MechanicsException(s.str());
            }
        }
    }
}

template<typename T>
int NuTo::BlockFullMatrix<T>::GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const
{
    int numCols = 0;
    for (auto dof : rDofTypes)
    {
        numCols += (*this)(dof, dof).GetNumColumns();
    }
    return numCols;
}



template<typename T>
int NuTo::BlockFullMatrix<T>::GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const
{
    int numRows = 0;
    for (auto dof : rDofTypes)
    {
        numRows += (*this)(dof, dof).GetNumRows();
    }
    return numRows;
}

template<typename T>
NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::BlockFullMatrix<T>::Export() const
{
    CheckDimensions();
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> result (GetNumActiveRows(), GetNumActiveColumns());

    int blockStartRow = 0;
    for (auto dofRow : activeDofTypes)
    {
        int blockStartCol = 0;
        for (auto dofCol : activeDofTypes)
        {
            const auto& subMatrix = (*this)(dofRow, dofCol);
            result.SetBlock(blockStartRow, blockStartCol, subMatrix);

            blockStartCol += subMatrix.GetNumColumns();
        }
        // every submatrix in the row dofRow has the same number of rows (CheckDimension)
        // --> one is picked and added to blockStartRow
        blockStartRow += (*this)(dofRow, dofRow).GetNumRows();
    }
    return result;
}

template<typename T>
NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::BlockFullMatrix<T>::Get(std::string rDofRow, std::string rDofCol) const
{
    return (*this)(Node::DofToEnum(rDofRow), Node::DofToEnum(rDofCol));
}

namespace NuTo
{
////! @brief stream operator for outputs with cout or files
template <typename T>
std::ostream& operator<< (std::ostream &rOut, const BlockFullMatrix<T>& rBlockMatrix)
{
    Eigen::IOFormat cleanFormat(Eigen::StreamPrecision, 0, " ", "\n", "|", " |");
    for (auto dof1 : rBlockMatrix.mDofStatus.GetActiveDofTypes())
    {
        for (auto dof2 : rBlockMatrix.mDofStatus.GetActiveDofTypes())
        {
            rOut << "[" << Node::DofToString(dof1) << " - " << Node::DofToString(dof2) << "]" << std::endl;
            rOut << rBlockMatrix(dof1,dof2).format(cleanFormat) << std::endl << std::endl;
        }
    }
    return rOut;
}
}//namespace NuTo

template std::ostream& NuTo::operator<< (std::ostream &rOut, const NuTo::BlockFullMatrix<double>& rBlockVector);



template class NuTo::BlockFullMatrix<double>;

#ifdef ENABLE_SERIALIZATION

namespace NuTo
{
template<> std::string BlockFullMatrix<double>::GetTypeId() const
{
    return "BlockFullMatrix<double>";
}
}

template void NuTo::BlockFullMatrix<double>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::BlockFullMatrix<double>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::BlockFullMatrix<double>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::BlockFullMatrix<double>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::BlockFullMatrix<double>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::BlockFullMatrix<double>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<typename T>
template<class Archive>
void NuTo::BlockFullMatrix<T>::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize BlockFullMarix" << "\n";
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(BlockStorageBase);
    ar& BOOST_SERIALIZATION_NVP(mData);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize BlockFullMarix \n";
#endif
}


BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::BlockFullMatrix<double>)
#endif //ENABLE_SERIALIZATION
