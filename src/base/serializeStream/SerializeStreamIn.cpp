#include "base/Exception.h"
#include "base/serializeStream/SerializeStreamIn.h"

NuTo::SerializeStreamIn::SerializeStreamIn(const std::string& rFile, bool rIsBinary)
    : SerializeStreamBase(rIsBinary)
{
    if (rIsBinary)
        mFileStream.open(rFile, std::ios_base::in | std::ios_base::binary);
    else
        mFileStream.open(rFile, std::ios_base::in);

    if (not mFileStream.is_open())
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Unable to read file " + rFile + ".");
}

void NuTo::SerializeStreamIn::Separator()
{
    if (mIsBinary)
    {
        char data;
        mFileStream.read(reinterpret_cast<char*>(&data), sizeof(SEPARATOR));
        if (data != *SEPARATOR)
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                  "Expected Separator \'" + std::string(SEPARATOR) + "\', got \'" + data + "\'");
    }
    else
    {
        std::string line;
        std::getline(mFileStream, line);
        if (line != std::string(SEPARATOR))
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                  "Expected Separator " + std::string(SEPARATOR) + ", got " + line);
    }
}
