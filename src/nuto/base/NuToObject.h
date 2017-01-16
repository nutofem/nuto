// $Id$

#pragma once

#include "nuto/base/Exception.h"

// STL:
#include <string>
#include <fstream>



#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //serialization

namespace NuTo {
class NuToObject {

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif// ENABLE_SERIALIZATION

public:
    NuToObject()
	{
		mVerboseLevel = 0;
#ifdef SHOW_TIME
		//! @brief ... show for each executed command the time required for the execution
        mShowTime = true;
#endif
	}

	virtual ~NuToObject()
	{}

    //! @brief ... sets the verbose level
    //! @param rVerboseLevel ... verbose level
	void SetVerboseLevel(const int rVerboseLevel)
	{
		mVerboseLevel=(unsigned short) rVerboseLevel;
	}

    //! @brief ... returns the verbose level
    //! @return verbose level
	int GetVerboseLevel()const
	{
		return mVerboseLevel;
	}



        //! @brief ... sets the showtime option
        //! @param rShowTime ... show time option
        void SetShowTime(bool rShowTime)
        {
#ifdef SHOW_TIME
            mShowTime=rShowTime;
#else
            (void)rShowTime;
#endif
        }

        //! @brief ... returns the show time optionl
        //! @return show time
        bool GetShowTime()const
        {
#ifdef SHOW_TIME
            return mShowTime;
#else
            return false;
#endif
        }

#ifdef SHOW_TIME
	//! @brief ...get the difference for exact time function in timespec format for testing purpose
	//! @param start, end ... time of begin and end
	//! @return  ... difference of times
	timespec diff(timespec start, timespec end)
	{
		timespec temp;
		if ((end.tv_nsec-start.tv_nsec)<0) {
			temp.tv_sec = end.tv_sec-start.tv_sec-1;
			temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
		} else {
			temp.tv_sec = end.tv_sec-start.tv_sec;
			temp.tv_nsec = end.tv_nsec-start.tv_nsec;
		}
		return temp;
	}


#endif

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
	virtual void Info()const=0;

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of NuToObject" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_NVP(mVerboseLevel);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of NuToObject" << std::endl;
#endif
    }

    //! @brief ... save the object to a file
	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
	virtual void Save (const std::string &filename, std::string rType )const
	{
	    try
	    {
	        //transform to uppercase
	        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
	        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
	        std::string tmpStr ( GetTypeId() );
	        std::string baseClassStr = tmpStr.substr ( 4,100 );
	        if (rType=="BINARY")
	        {
	            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
	            oba & boost::serialization::make_nvp ( "Object_type", tmpStr );
	            oba & boost::serialization::make_nvp(tmpStr.c_str(), *this);
	        }
	        else if (rType=="XML")
	        {
	            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
	            oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
	            oxa & boost::serialization::make_nvp(tmpStr.c_str(), *this);
	        }
	        else if (rType=="TEXT")
	        {
	            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
	            ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
	            ota & boost::serialization::make_nvp(tmpStr.c_str(), *this);
	        }
	        else
	        {
	            throw Exception ( "[NewmarkDirect::Save]File type not implemented." );
	        }
	    }
	    catch ( boost::archive::archive_exception e )
	    {
	        std::string s ( std::string ( "[NewmarkDirect::Save]File save exception in boost - " ) +std::string ( e.what() ) );
	        std::cout << s << "\n";
	        throw Exception ( s );
	    }
	    catch ( Exception &e )
	    {
	        throw e;
	    }
	    catch ( std::exception &e )
	    {
	        throw Exception ( e.what() );
	    }
	    catch ( ... )
	    {
	        throw Exception ( "[NewmarkDirect::Save]Unhandled exception." );
	    }
	}

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    virtual void Restore (const std::string &filename, std::string rType )
    {
        try
        {
            //transform to uppercase
            std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
            std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
            std::string tmpString;
            if (rType=="BINARY")
            {
                boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
                oba & boost::serialization::make_nvp ( "Object_type", tmpString );
                if ( tmpString!=GetTypeId() )
                    throw Exception ( "[NewmarkDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
                oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
            }
            else if (rType=="XML")
            {
                boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
                oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
                if ( tmpString!=GetTypeId() )
                    throw Exception ( "[NewmarkDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
                oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
            }
            else if (rType=="TEXT")
            {
                boost::archive::text_iarchive ota ( ifs, std::ios::binary );
                ota & boost::serialization::make_nvp ( "Object_type", tmpString );
                if ( tmpString!=GetTypeId() )
                    throw Exception ( "[NewmarkDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
                ota & boost::serialization::make_nvp(tmpString.c_str(), *this);
            }
            else
            {
                throw Exception ( "[Matrix::Restore]File type not implemented" );
            }
        }
        catch ( Exception &e )
        {
            throw e;
        }
        catch ( std::exception &e )
        {
            throw Exception ( e.what() );
        }
        catch ( ... )
        {
            throw Exception ( "[NewmarkDirect::Restore]Unhandled exception." );
        }
    }

#endif //serialization

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId()const=0;

protected:
    unsigned short mVerboseLevel; //!< verbose level between 0 (no output) and 10 (all actions are commented on the console)

#ifdef SHOW_TIME
    //! @brief ... show for each executed command the time required for the execution
    bool mShowTime;
#endif
};
}