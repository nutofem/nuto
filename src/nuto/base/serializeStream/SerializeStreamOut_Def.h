#pragma once

#include "nuto/base/serializeStream/SerializeStreamBase.h"

namespace Eigen
{
template <typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols> class Eigen::Matrix;
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

    //! @brief perfect(?) forwarding to WriteData
    template<class ...Args>
    void SerializeData(const Args&... rValue);

    //! @brief Writes a double to the stream
    void WriteData(double rValue);

    //! @brief Writes an Eigen::Matrix to the stream
    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void WriteData(const Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols> & rMatrix);
};
}
