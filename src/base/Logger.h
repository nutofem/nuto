#pragma once

#include <iostream>
#include <fstream>
#include <string>

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date July 2011
//! @brief ... logger class for redirecting output to different locations/files
class Logger
{

public:
    //! @brief ...constructor
    Logger();

    //! @brief ..opens the file stored in the string mLogFileName
    void OpenFile();

    //! @brief ..open file
    //! @param rFileString file string to open
    void OpenFile(std::string rFileString);

    //! @brief ..Close file
    void CloseFile();

    //! @brief ..sets the output to be forwarded to the console (false) or only to the file (true)
    //! @param rQuiet see explanation above
    void SetQuiet(bool rQuiet);

    //! @brief ..Writes a message to the log and to console
    template <typename T>
    void Out(const T& rObject, bool rEof)
    {
        if (!mQuiet)
        {
            std::cout << rObject << (rEof ? "\n" : "");
            std::cout.flush();
        }
        if (mLogFile.is_open())
        {
            mLogFile << rObject << (rEof ? "\n" : "");
            mLogFile.flush();
        }
    }

    //! @brief Info routine that prints general information about the object (detail according to verbose level)
    void Info() const
    {
        std::cout << "LogFileName " << mLogFileName << " is quiet " << mQuiet << std::endl;
    }


protected:
    std::ofstream mLogFile; //!< Logfile for output.
    std::string mLogFileName ; //!< LogfileName for output.
    bool mQuiet = false; //!< If true, no writing to console.;
};

//! Generic output command.
template <typename T>
Logger& operator<<(Logger& rLogger, const T& rObject)
{
    rLogger.Out(rObject, false);
    return rLogger;
}
} // namespace NuTo
