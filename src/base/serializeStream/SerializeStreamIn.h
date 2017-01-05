#pragma once

#include "base/serializeStream/SerializeStreamBase.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{
//! @brief Serialize input stream
class SerializeStreamIn : public SerializeStreamBase
{
public:

    //! @brief ctor
    //! @param rFile file name
    //! @param rIsBinary flag to enable/disable binary file read
    SerializeStreamIn(const std::string& rFile, bool rIsBinary);

    virtual ~SerializeStreamIn()                     = default;
    SerializeStreamIn(const SerializeStreamIn&)      = delete;

    //! @brief in-stream operator
    //! @param rStream NuTo input stream
    //! @param rData the data from the input stream goes here
    template <typename T>
    inline friend SerializeStreamIn& operator >> (SerializeStreamIn& rStream, T& rData)
    {
        rStream.Serialize(rData);
        return rStream;
    }

    template <typename T>
    void Serialize(T& rData)        {rData.NuToSerializeLoad(*this);}
    void Serialize(double& rData)   {SerializePrimitiveType(rData);}
    void Serialize(int& rData)      {SerializePrimitiveType(rData);}
    void Serialize(bool& rData)     {SerializePrimitiveType(rData);}

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void Serialize(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix)
    {
            LoadMatrix(rMatrix);
    }

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void LoadMatrix(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix)
    {
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

    //! @brief reads a sequence from the stream
    void Separator();

private:
    template <typename T>
    void SerializePrimitiveType(T& rData)
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

};
} // namespace NuTo
