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



namespace NuTo
{

template <typename T>
SerializeStreamOut& operator<<(SerializeStreamOut& rStream, const T& rData)
{
    rData.NuToSerializeSave(rStream);
    return rStream;
}

SerializeStreamOut& operator<<(SerializeStreamOut& rStream, double rData)
{
    if (rStream.mIsBinary)
    {
        rStream.mFileStream.write(reinterpret_cast<const char*>(&rData), sizeof(double));
    }
    else
    {
        rStream.mFileStream << typeid(rData).name() << ":" << std::endl; // one line of debug info
        rStream.mFileStream << rData << std::endl;
    }
    return rStream;
}

template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
SerializeStreamOut& operator<<(SerializeStreamOut& rStream, const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix)
{
    rStream.SaveMatrix(rMatrix);
    return rStream;
}

}   // namespace NuTo

template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
void NuTo::SerializeStreamOut::SaveMatrix(const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix)
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
};
