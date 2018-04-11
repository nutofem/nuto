#pragma once

#include <iostream>
#include <fstream>
#include <string>

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date July 2011
//! @brief logger class for redirecting output to different locations/files
class Logger
{

public:
    //! @brief ...constructor
    Logger(std::string prefix = "", bool isQuiet = true);

    //! @brief ..open file
    //! @param filename file string to open
    void OpenFile(std::string filename);

    //! @brief ..Close file
    void CloseFile();

    //! @brief ..sets the output to be forwarded to the console (false) or only to the file (true)
    //! @param isQuiet see explanation above
    void SetQuiet(bool isQuiet);

    //! @brief ..Writes a message to the log and to console
    template <typename T>
    void Out(const T& rObject, bool withNewline = false)
    {
        std::string optionalPrefix = OptionalPrefix();
        if (!mIsQuiet)
        {
            std::cout << optionalPrefix << rObject;
            std::cout.flush();
        }
        if (mLogFile.is_open())
        {
            mLogFile << optionalPrefix << rObject;
            mLogFile.flush();
        }
        mPrintPrefix = withNewline;
    }

    std::string OptionalPrefix();

    //! @brief Info routine that prints general information about the object (detail according to verbose level)
    void Info() const
    {
        std::cout << "LogFileName " << mLogFileName << " is quiet " << mIsQuiet << std::endl;
    }


protected:
    std::ofstream mLogFile; //!< Logfile for output.
    std::string mLogFileName; //!< LogfileName for output.
    bool mIsQuiet = false; //!< If true, no writing to console.;

    std::string mPrefix; //!< prefix written on the beginning of the log line
    bool mPrintPrefix = true;
};

//! Generic output command.
template <typename T>
inline Logger& operator<<(Logger& rLogger, const T& t)
{
    rLogger.Out(t);
    return rLogger;
}

Logger& operator<<(Logger& rLogger, const char& t);
Logger& operator<<(Logger& rLogger, const std::string& t);
Logger& operator<<(Logger& rLogger, const char* t);

//! @brief set of predefined global loggers
struct Log
{
    static Logger Debug; // prefix "Debug| "
    static Logger Info; // prefix "Info|  "
    static Logger Error; // prefix "Error| "
};
} // namespace NuTo
