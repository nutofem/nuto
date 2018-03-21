#include "nuto/base/Logger.h"
#include "nuto/base/Exception.h"

NuTo::Logger::Logger()
    : mLogFile()
    , mLogFileName()
{
}

void NuTo::Logger::OpenFile()
{
    mLogFile.close();
    mLogFile.clear();
    if (mLogFileName.length() != 0)
    {
        mLogFile.open(mLogFileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
        if (!mLogFile.is_open())
            throw Exception("[NuTo::Logger::OpenFile]File could not be opened.");
    }
}

void NuTo::Logger::OpenFile(std::string rFileString)
{
    mLogFileName = rFileString;
    mLogFile.close();
    mLogFile.clear();
    mLogFile.open(rFileString.c_str());
    if (!mLogFile.is_open())
        throw Exception("[NuTo::Logger::OpenFile] File " + std::string(rFileString.c_str()) + " could not be opened.");
}

void NuTo::Logger::CloseFile()
{
    mLogFile.close();
}

void NuTo::Logger::SetQuiet(bool rQuiet)
{
    mQuiet = rQuiet;
}
