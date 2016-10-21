#include "nuto/base/serializeStream/SerializeStreamIn.h"

NuTo::SerializeStreamIn::SerializeStreamIn(const std::string& rFile, bool rIsBinary)
    : SerializeStreamBase(rIsBinary)
{
    if (rIsBinary)
        mFileStream.open(rFile, std::ios_base::in | std::ios_base::binary);
    else
        mFileStream.open(rFile, std::ios_base::in);
}