#pragma once
#include <fstream>

namespace NuTo
{
//! @brief Base class for the NuTo SerializeStream
//! Used for the serialization of data values only. The class structure serialization is done with boost::serialize
//! @remark the std::fstream is added as a member function (not derived) to avoid unintended behavior with its
//! operator<<
class SerializeStreamBase
{
public:
    SerializeStreamBase(bool rIsBinary);
    virtual ~SerializeStreamBase() = default;
    SerializeStreamBase(const SerializeStreamBase&) = delete;

protected:
    const bool mIsBinary;
    std::fstream mFileStream;
    static constexpr const char* SEPARATOR = "#";
};
}
