#pragma once

#include "nuto/base/serializeStream/SerializeStreamIn_Def.h"
#include <eigen3/Eigen/Core>

NuTo::SerializeStreamIn::SerializeStreamIn(const std::string& rFile, bool rIsBinary)
    : SerializeStreamBase(rIsBinary)
{
    if (rIsBinary)
        mFileStream.open(rFile, std::ios_base::in | std::ios_base::binary);
    else
        mFileStream.open(rFile, std::ios_base::in);
}

template <typename T>
void NuTo::SerializeStreamIn::NuToSerializeNumber(T &rData)
{
    if (mIsBinary)
    {
        mFileStream.read(reinterpret_cast<char*>(&rData), sizeof(T));
    }
    else
    {
        std::string line;
        std::getline(mFileStream, line); // ignore one line of debug info
        std::getline(mFileStream, line); // extract value
        rData = static_cast<T>(std::stod(line));
    }
}


template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
void NuTo::SerializeStreamIn::NuToSerializeMatrix(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix) {
    const auto &rows = rMatrix.rows();
    const auto &cols = rMatrix.cols();
    auto *data = rMatrix.data();
    if (mIsBinary) {
        mFileStream.read(reinterpret_cast<char *>(data), rows * cols * sizeof(T));
    } else {
        std::string line;
        std::getline(mFileStream, line); // ignore one line of debug info
        for (int i = 0; i < rows * cols; ++i) {
            std::getline(mFileStream, line);
            data[i] = std::stod(line);
        }
    }
}