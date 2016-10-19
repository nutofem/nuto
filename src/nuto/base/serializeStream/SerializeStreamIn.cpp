#include "nuto/base/serializeStream/SerializeStreamIn.h"

NuTo::SerializeStreamIn::SerializeStreamIn(const std::string& rFile, bool rIsBinary)
    : SerializeStreamBase(rIsBinary)
{
    if (rIsBinary)
        mFileStream.open(rFile, std::ios_base::in | std::ios_base::binary);
    else
        mFileStream.open(rFile, std::ios_base::in);
}

void NuTo::SerializeStreamIn::Serialize(double& rData)
{
    if (mIsBinary)
    {
        mFileStream.read(reinterpret_cast<char*>(&rData), sizeof(double));
    }
    else
    {
        std::string line;
        std::getline(mFileStream, line); // ignore one line of debug info
        std::getline(mFileStream, line); // extract value
        rData = std::stod(line);
    }
}


void NuTo::SerializeStreamIn::Serialize(bool& rData)
{
    if (mIsBinary)
    {
        mFileStream.read(reinterpret_cast<char*>(&rData), sizeof(bool));
    }
    else
    {
        std::string line;
        std::getline(mFileStream, line); // ignore one line of debug info
        std::getline(mFileStream, line); // extract value
        rData = static_cast<bool>(std::stoi(line));
    }
}
