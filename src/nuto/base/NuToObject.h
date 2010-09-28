// $Id: $

#ifndef NuToObject_H
#define NuToObject_H

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //serialization

namespace NuTo {
class NuToObject {
public:
	NuToObject()
	{
		mVerboseLevel = 0;
#ifdef SHOW_TIME
    //! @brief ... show for each executed command the time required for the execution
        mShowTime = false;
#endif
		mVerboseLevel = 0;
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

#ifdef SHOW_TIME
	//! @brief ... sets the showtime option
    //! @param rShowTime ... show time option
	void SetShowTime(bool rShowTime)
	{
		mShowTime=rShowTime;
	}

    //! @brief ... returns the show time optionl
    //! @return show time
	bool GetShowTime()const
	{
		return mShowTime;
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
	virtual void Save (const std::string &filename, std::string rType )const=0;
/*	{
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
				///    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Matrix<T>)
				// work around for BOOST_SERIALIZATION_BASE_OBJECT_NVP
				oba & boost::serialization::make_nvp ( baseClassStr.c_str(),boost::serialization::base_object< Matrix<T> > ( *this ) );
				std::vector<T> dataVec ( GetNumRows() *GetNumColumns() );
				int numRows = mEigenMatrix.rows(),
						numColumns = mEigenMatrix.cols();
				memcpy ( & ( dataVec[0] ),mEigenMatrix.data(),GetNumRows() *GetNumColumns() *sizeof ( T ) );
				oba & BOOST_SERIALIZATION_NVP ( dataVec )
				& BOOST_SERIALIZATION_NVP ( numRows )
				& BOOST_SERIALIZATION_NVP ( numColumns );
			}
			else if (rType=="XML")
			{
				boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
				oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
				// work around for BOOST_SERIALIZATION_BASE_OBJECT_NVP
				//oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Matrix<T>);
				oxa & boost::serialization::make_nvp ( baseClassStr.c_str(),boost::serialization::base_object< Matrix<T> > ( *this ) );
				std::vector<T> dataVec ( GetNumRows() *GetNumColumns() );
				int numRows = mEigenMatrix.rows(),
						numColumns = mEigenMatrix.cols();
				memcpy ( & ( dataVec[0] ),mEigenMatrix.data(),GetNumRows() *GetNumColumns() *sizeof ( T ) );
				oxa & BOOST_SERIALIZATION_NVP ( dataVec )
				& BOOST_SERIALIZATION_NVP ( numRows )
				& BOOST_SERIALIZATION_NVP ( numColumns );
			}
			else if (rType=="TEXT")
			{
				boost::archive::text_oarchive ota ( ofs, std::ios::binary );
				ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
				// work around for BOOST_SERIALIZATION_BASE_OBJECT_NVP
				//ota & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Matrix<T>);
				ota & boost::serialization::make_nvp ( baseClassStr.c_str(),boost::serialization::base_object< Matrix<T> > ( *this ) );
				std::vector<T> dataVec ( GetNumRows() *GetNumColumns() );
				int numRows = mEigenMatrix.rows(),
						numColumns = mEigenMatrix.cols();
				memcpy ( & ( dataVec[0] ),mEigenMatrix.data(),GetNumRows() *GetNumColumns() *sizeof ( T ) );
				ota & BOOST_SERIALIZATION_NVP ( dataVec )
				& BOOST_SERIALIZATION_NVP ( numRows )
				& BOOST_SERIALIZATION_NVP ( numColumns );
			}
			else
			{
				throw MathException ( "[FullMatrix::Save]File type not implemented." );
			}
		}
		catch ( boost::archive::archive_exception e )
		{
			std::string s ( std::string ( "[FullMatrix::Save]File save exception in boost - " ) +std::string ( e.what() ) );
			std::cout << s << "\n";
			throw MathException ( s );
		}
		catch ( MathException &e )
		{
			throw e;
		}
		catch ( std::exception &e )
		{
			throw MathException ( e.what() );
		}
		catch ( ... )
		{
			throw MathException ( "[Matrix::Save]Unhandled exception." );
		}
    }
    */

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    virtual void Restore (const std::string &filename, std::string rType )=0;
/*
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
					throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

				// work around for BOOST_SERIALIZATION_BASE_OBJECT_NVP
				//    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Matrix<T>)
				std::string baseClassStr = tmpString.substr ( 4,100 );
				oba & boost::serialization::make_nvp ( baseClassStr.c_str(),boost::serialization::base_object< Matrix<T> > ( *this ) );

				std::vector<T> dataVec;
				int numRows, numColumns;
				oba & BOOST_SERIALIZATION_NVP ( dataVec )
				& BOOST_SERIALIZATION_NVP ( numRows )
				& BOOST_SERIALIZATION_NVP ( numColumns );
				mEigenMatrix.resize ( numRows,numColumns );
				memcpy ( mEigenMatrix.data(),& ( dataVec[0] ),numRows*numColumns*sizeof ( T ) );
			}
			else if (rType=="XML")
			{
				boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
				oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
				if ( tmpString!=GetTypeId() )
					throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

				// work around for BOOST_SERIALIZATION_BASE_OBJECT_NVP
				//    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Matrix<T>)
				std::string baseClassStr = tmpString.substr ( 4,100 );
				oxa & boost::serialization::make_nvp ( baseClassStr.c_str(),boost::serialization::base_object< Matrix<T> > ( *this ) );

				std::vector<T> dataVec;
				int numRows, numColumns;
				oxa & BOOST_SERIALIZATION_NVP ( dataVec )
				& BOOST_SERIALIZATION_NVP ( numRows )
				& BOOST_SERIALIZATION_NVP ( numColumns );
				mEigenMatrix.resize ( numRows,numColumns );
				memcpy ( mEigenMatrix.data(),& ( dataVec[0] ),numRows*numColumns*sizeof ( T ) );
			}
			else if (rType=="TEXT")
			{
				boost::archive::text_iarchive ota ( ifs, std::ios::binary );
				ota & boost::serialization::make_nvp ( "Object_type", tmpString );
				if ( tmpString!=GetTypeId() )
					throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

				// work around for BOOST_SERIALIZATION_BASE_OBJECT_NVP
				//    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Matrix<T>)
				std::string baseClassStr = tmpString.substr ( 4,100 );
				ota & boost::serialization::make_nvp ( baseClassStr.c_str(),boost::serialization::base_object< Matrix<T> > ( *this ) );

				std::vector<T> dataVec;
				int numRows, numColumns;
				ota & BOOST_SERIALIZATION_NVP ( dataVec )
				& BOOST_SERIALIZATION_NVP ( numRows )
				& BOOST_SERIALIZATION_NVP ( numColumns );
				mEigenMatrix.resize ( numRows,numColumns );
				memcpy ( mEigenMatrix.data(),& ( dataVec[0] ),numRows*numColumns*sizeof ( T ) );
			}
			else
			{
				throw MathException ( "[Matrix::Restore]File type not implemented" );
			}
		}
		catch ( MathException &e )
		{
			throw e;
		}
		catch ( std::exception &e )
		{
			throw MathException ( e.what() );
		}
		catch ( ... )
		{
			throw MathException ( "[Matrix::Restore]Unhandled exception." );
		}
    }
*/
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
#endif //NuToObject_H
