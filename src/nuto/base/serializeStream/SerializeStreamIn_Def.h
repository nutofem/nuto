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

    //! @brief perfect(?) forwarding to ReadData
    template<class ...Args>
    void SerializeData(Args&... rValue);

    //! @brief Reads a double from the stream
    void ReadData(double& rValue);

    //! @brief Reads an Eigen::Matrix from the stream
    template<typename T, int TRows, int TCols, int TOptions, int TMaxRows, int TMaxCols>
    void ReadData(Eigen::Matrix<T, TRows, TCols, TOptions, TMaxRows, TMaxCols> & rMatrix);
};
}
