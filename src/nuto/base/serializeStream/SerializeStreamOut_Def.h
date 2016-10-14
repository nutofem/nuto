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

    SerializeStreamOut(const SerializeStreamOut&)     = delete;

    template <typename T>
    friend SerializeStreamOut& operator<<(SerializeStreamOut& rStream, const T& rData)
    {
        rData.NuToSerializeWrite(rStream);
        return rStream;
    }

    template <typename T>
    void NuToSerializeNumber(const T &rData);

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void NuToSerializeMatrix(const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);
};
} // namespace NuTo
