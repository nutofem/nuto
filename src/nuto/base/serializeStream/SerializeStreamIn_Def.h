#pragma once

#include "nuto/base/serializeStream/SerializeStreamBase.h"

namespace Eigen
{
template <typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols> class Matrix;
}

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
    void Serialize(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void LoadMatrix(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);

private:
    template <typename T>
    void SerializePrimitiveType(T& rData);

};
} // namespace NuTo
