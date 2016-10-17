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
    friend SerializeStreamOut& operator<<(SerializeStreamOut& rStream, const T& rData);

    friend SerializeStreamOut& operator<<(SerializeStreamOut& rStream, double rData);

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    friend SerializeStreamOut& operator<<(SerializeStreamOut& rStream, const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);

    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void SaveMatrix(const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols>& rMatrix);

};
} // namespace NuTo
