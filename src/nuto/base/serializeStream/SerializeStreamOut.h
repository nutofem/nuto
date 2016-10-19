#pragma once

#include "nuto/base/serializeStream/SerializeStreamOut_Def.h"
#include <eigen3/Eigen/Core>


template <typename T>
void NuTo::SerializeStreamOut::Serialize(T& rData)
{
    rData.NuToSerializeSave(*this);
}

template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
void NuTo::SerializeStreamOut::Serialize(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix)
{
    SaveMatrix(rMatrix);
}


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
}
