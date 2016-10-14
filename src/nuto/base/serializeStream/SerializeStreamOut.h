#pragma once

#include "nuto/base/serializeStream/SerializeStreamOut_Def.h"
#include <eigen3/Eigen/Core>
#include <typeinfo>

NuTo::SerializeStreamOut::SerializeStreamOut(const std::string& rFile, bool rIsBinary)
    : SerializeStreamBase(rIsBinary)
{
    if (rIsBinary)
    {
        mFileStream.open(rFile, std::ios_base::out | std::ios_base::binary);
    }
    else
    {
        constexpr int doublePrecision = 17; // somehow required precision
        mFileStream.open(rFile, std::ios_base::out);
        mFileStream.precision(doublePrecision);
    }
}

template <typename T>
void NuTo::SerializeStreamOut::NuToSerializeNumber(const T &rData)
{
    if (mIsBinary)
    {
        mFileStream.write(reinterpret_cast<const char*>(&rData), sizeof(T));
    }
    else
    {
        mFileStream << typeid(rData).name() << ":" << std::endl; // one line of debug info
        mFileStream << rData << std::endl;
    }
}

template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
void NuTo::SerializeStreamOut::NuToSerializeMatrix(const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols> &rMatrix)
{
    const auto& rows = rMatrix.rows();
    const auto& cols = rMatrix.cols();
    const auto& data = rMatrix.data();
    if (mIsBinary)
    {
        mFileStream.write(reinterpret_cast<const char*>(data), rows*cols*sizeof(T));
    }
    else
    {
        // one line of debug info:
        mFileStream << "Matrix ( " << rows << " x " << cols << " ): " << std::endl;
        for (int i = 0; i < rows*cols; ++i)
            mFileStream << data[i] << std::endl;
    }
}