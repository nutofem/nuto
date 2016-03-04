#ifndef NUTO_LOGGER_H
#define NUTO_LOGGER_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/export.hpp>
#include <boost/serialization/access.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/NuToObject.h"
#include "nuto/math/FullMatrix_Def.h"

#include <iostream>
#include <fstream>
#include <string>

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date July 2011
//! @brief ... logger class for redirecting output to different locations/files
class Logger;
class Logger : public NuToObject
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief ...constructor
    Logger();

#ifdef ENABLE_SERIALIZATION
// serializes the class
template<class Archive>
void load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of Logger" << std::endl;
#endif
    bool wasOpen = false;
    ar & BOOST_SERIALIZATION_NVP(mQuiet);
    ar & BOOST_SERIALIZATION_NVP(mLogFileName);
    ar & boost::serialization::make_nvp("isOpen", wasOpen);
    if (wasOpen)
    {
        //make sure it has been closed before, otherwise do nothing
        mLogFile.close();
        mLogFile.clear();
        mLogFile.open(mLogFileName.c_str(),std::ios_base::app);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of Logger" << std::endl;
#endif
}

template<class Archive>
void save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of Logger" << std::endl;
#endif
    bool isOpen(mLogFile.is_open());
    ar & BOOST_SERIALIZATION_NVP(mQuiet);
    ar & BOOST_SERIALIZATION_NVP(mLogFileName);
    ar & boost::serialization::make_nvp("isOpen", isOpen);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of Logger" << std::endl;
#endif
}

BOOST_SERIALIZATION_SPLIT_MEMBER()

#endif //ENABLE_SERIALIZATION

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
    template<typename T> void Out(const T & rObject, bool rEof)
    {
        if (!mQuiet)
        {
            std::cout << rObject << (rEof ? "\n" : "");
            std::cout.flush();
        }
        if(mLogFile.is_open())
        {
        	mLogFile << rObject << (rEof ? "\n" : "");
            mLogFile.flush();
        }
    }

    //! @brief ..logs a NuTo::FullMatrix
    //! @param rInt1 parameters total number of digits to be plotted
    //! @param rInt2 parameters number of digits after the comma to be plotted
    void Out(const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rObject, int rInt1, int rInt2, bool rScientific=false);

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const
    {
        return std::string("Logger");
    }

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
	void Info()const
	{
		std::cout << "LogFileName " << mLogFileName << " is quiet " << mQuiet << std::endl;
	}


#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
	void Save (const std::string &filename, std::string rType )const;

	//! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore (const std::string &filename, std::string rType );
#endif //ENABLE_SERIALIZATION

protected:
    std::ofstream mLogFile;     //!< Logfile for output.
    std::string mLogFileName;   //!< LogfileName for output.
    int mQuiet;                 //!< If true, no writing to console.;
};

//! Generic output command.
template <typename T> Logger& operator<<(Logger &rLogger, const T &rObject)
{
	rLogger.Out(rObject, false);
    return rLogger;
}
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::Logger)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif  // NUTO_LOGGER_H
