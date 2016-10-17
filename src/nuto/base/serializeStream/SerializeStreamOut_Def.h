#pragma once

#include "nuto/base/serializeStream/SerializeStreamBase.h"

namespace Eigen
{
template <typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols> class Matrix;
}

namespace NuTo
{
//! @brief Serialize output stream
class SerializeStreamOut : public SerializeStreamBase
{
public:

    //! @brief ctor
    //! @param rFile file name
    //! @param rIsBinary flag to enable/disable binary file read
    SerializeStreamOut(const std::string& rFile, bool rIsBinary);
    virtual ~SerializeStreamOut()                     = default;

    template <typename T>
    inline friend SerializeStreamOut& operator << (SerializeStreamOut& rStream, T& rData)
    {
        rStream.Serialize(rData);
        return rStream;
    }

    template <typename T>
    void Serialize(T& rData);

    void Serialize(double rData);

    void Serialize(bool rData);

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void Serialize(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void SaveMatrix(const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);

};
} // namespace NuTo
