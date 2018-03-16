#include "nuto/base/serializeStream/SerializeStreamOut.h"
#include "nuto/base/Exception.h"

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
    if (not mFileStream.is_open())
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Unable to create file " + rFile + ".");
}

void NuTo::SerializeStreamOut::Separator()
{
    if (mIsBinary)
    {
        mFileStream.write(SEPARATOR, sizeof(SEPARATOR));
    }
    else
    {
        mFileStream << SEPARATOR << '\n';
    }
}
