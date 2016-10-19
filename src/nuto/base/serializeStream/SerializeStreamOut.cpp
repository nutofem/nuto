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


void NuTo::SerializeStreamOut::Serialize(double rData)
{
    if (mIsBinary)
    {
        mFileStream.write(reinterpret_cast<const char*>(&rData), sizeof(double));
    }
    else
    {
        mFileStream << "double:" << std::endl; // one line of debug info
        mFileStream << rData << std::endl;
    }
}

void NuTo::SerializeStreamOut::Serialize(bool rData)
{
    if (mIsBinary)
    {
        mFileStream.write(reinterpret_cast<const char*>(&rData), sizeof(bool));
    }
    else
    {
        mFileStream << "bool:" << std::endl; // one line of debug info
        mFileStream << static_cast<int>(rData) << std::endl; // make sure its 0 or 1
    }
}


