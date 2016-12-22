#include "base/Logger.h"
#include "base/Exception.h"

NuTo::Logger::Logger()
{
	mQuiet = false;
}

#ifdef ENABLE_SERIALIZATION
//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Logger::Save (const std::string &filename, std::string rType )const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw Exception("[NuTo::Logger::Save] Error opening file.");
        }

        // write data to file
        std::string typeIdString(this->GetTypeId());
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            oba & boost::serialization::make_nvp ("Object_type", typeIdString );
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            std::string tmpString(this->GetTypeId());
            oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp("Object_type", typeIdString );
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw Exception ( "[NuTo::Logger::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception& e )
    {
        std::string s ( std::string ( "[NuTo::Logger::Save]File save exception in boost - " ) + std::string ( e.what() ) );
        throw Exception ( s );
    }
    catch (Exception &e )
    {
        throw;
    }
    catch ( std::exception &e )
    {
        throw Exception ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Logger::Save] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Logger::Restore (const std::string &filename, std::string rType )
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw Exception("[NuTo::Logger::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw Exception ( "[NuTo::Logger::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw Exception ( "[NuTo::Logger::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw Exception ( "[NuTo::Logger::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw Exception ( "[NuTo::Logger::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception& e )
    {
        std::string s ( std::string ( "[NuTo::Logger::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
        throw Exception ( s );
    }
    catch ( Exception &e )
    {
        throw;
    }
    catch ( std::exception &e )
    {
        throw Exception ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Logger::Restore] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

#endif //SERIALIZATION

//! @brief ..opens the file stored in the string mLogFileName
void NuTo::Logger::OpenFile()
{
    mLogFile.close();
    mLogFile.clear();
    if (mLogFileName.length()!=0)
    {
		mLogFile.open(mLogFileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		if (!mLogFile.is_open())
			throw Exception ( "[NuTo::Logger::OpenFile]File could not be opened." );
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
		throw Exception ( "[NuTo::Logger::OpenFile] File " + std::string(rFileString.c_str()) + " could not be opened." );
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
    mQuiet=rQuiet;
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Logger)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
