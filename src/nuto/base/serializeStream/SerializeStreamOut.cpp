#include "nuto/base/serializeStream/SerializeStreamOut.h"

NuTo::SerializeStreamOut::SerializeStreamOut(const std::string& rFile, bool rIsBinary)
    : SerializeStreamBase(rIsBinary)
{
    if (rIsBinary)
    {
        mFileStream.open(rFile, std::ios_base::out | std::ios_base::binary);
    }
    else
    {
        constexpr int doublePrecision = 17; // somehow required precision
        mFileStream.open(rFile, std::ios_base::out);
        mFileStream.precision(doublePrecision);
    }
}