#include "nuto/base/Logger.h"
#include "nuto/base/Exception.h"

NuTo::Logger::Logger()
    : mLogFile()
    , mLogFileName()
{
}

//! @brief ..opens the file stored in the string mLogFileName
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

//! @brief ..open file
//! @param rFileString file string to open
void NuTo::Logger::OpenFile(std::string rFileString)
{
    mLogFileName = rFileString;
    mLogFile.close();
    mLogFile.clear();
    mLogFile.open(rFileString.c_str());
    if (!mLogFile.is_open())
        throw Exception("[NuTo::Logger::OpenFile] File " + std::string(rFileString.c_str()) + " could not be opened.");
}

//! @brief ..reset file and open it again
//! @param rFileString file string to open
void NuTo::Logger::CloseFile()
{
    mLogFile.close();
}

//! @brief ..sets the output to be forwarded to the console (false) or only to the file (true)
//! @param rQuiet see explanation above
void NuTo::Logger::SetQuiet(bool rQuiet)
{
    mQuiet = rQuiet;
}
