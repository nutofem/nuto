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

    template <typename T>
    void NuToSerializeNumber(T &rData);

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void NuToSerializeMatrix(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);


    //! @brief in-stream operator
    //! @param rStream NuTo input stream
    //! @param rData the data from the input stream goes here
    template <typename T>
    friend SerializeStreamIn& operator>>(SerializeStreamIn& rStream, T& rData)
    {
        rData.NuToSerializeRead(rStream);
        return rStream;
    }
};
} // namespace NuTo
