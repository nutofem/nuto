// $Id$

#ifndef FULL_MATRIX_H
#define FULL_MATRIX_H

//plugins for the eigen matrices
#define EIGEN_MATRIXBASE_PLUGIN "MatrixBaseAddons.h"

#include <iomanip>  //for setw
#include <fstream>  //for file acces
#include <vector>
#include <string>
#include <iostream>
//using namespace std;
#include <boost/tokenizer.hpp>


#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#endif  // ENABLE_SERIALIZATION

#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>

#include "nuto/math/FullMatrix_Def.h"
#include "nuto/math/FullVector.h"

#include "nuto/base/Logger.h"
#include "nuto/math/NuToMath.h"
#include "nuto/math/MathException.h"
#include "nuto/math/fortran_routines.h"
#include "nuto/math/SparseMatrix.h"


namespace NuTo
{
//! @brief ... default constructor
template <class T, int rows, int cols>
FullMatrix<T,rows,cols>::FullMatrix()
{
}

//! @brief ... constructor
//! @param rNumRows ... number of columns
//! @param rNumColumns ... number of rows
template <class T, int rows, int cols>
FullMatrix<T,rows,cols>::FullMatrix( int rNumRows, int rNumColumns)
{
	if ( rNumRows*rNumColumns>0 )
		this->setZero ( rNumRows,rNumColumns );
	else
		this->resize ( rNumRows,rNumColumns );
}

// convert numpy array python to eign matrix c++
#ifdef ENABLE_NUMPY
	template <class T, int rows, int cols>
	void FullMatrix<T,rows,cols>::convrtNumpyToMatrix(T * inData,int rRow, int rCol)

	{
		int count =0;
		for(int i=0;i< rRow;i++)
		{
		   for(int j=0;j<rCol;j++)
			{
			(*this)(i,j)=inData[count];
			count++;
			}
		}
	}

	// convert eign matrix from c++ to numpy array python
	template <class T, int rows, int cols>
	void FullMatrix<T,rows,cols>::convrtMatrixToNumpy(T * indata,int rRow,int rCol)
	{
		int count=0;
		for(int i=0;i<rRow;i++)
		{
			for(int j=0;j<rCol;j++)
			{
			indata[count]=(*this)(i,j) ;
			count++;
			}
		}
	}
#endif //ENABLE_NUMPY


//! @brief ... constructor
//! @param rEigenMatrix ... other matrix
template<class T, int rows, int cols>
template<typename OtherDerived>
FullMatrix<T,rows,cols>::FullMatrix ( const Eigen::MatrixBase<OtherDerived>& rEigenMatrix): Eigen::Matrix<T,rows,cols>(rEigenMatrix)
{
}

    //! @brief ... constructur
    //! @param numRows_ ... number of columns
    //! @param numColumns_ ... number of rows
    //! @param entries_ ... vector containing the matrix in column-major orientation
template<class T, int rows, int cols>
FullMatrix<T,rows,cols>::FullMatrix ( int numRows_, int numColumns_, const std::vector<T>& entries_ )
{
	this->resize ( numRows_,numColumns_ );
	const T *ptr = &entries_[0];
	for ( int j=0; j<this->cols(); ++j )              // loop over columns
		for ( int i=0; i<this->rows(); ++i, ptr++ )   // loop over rows
			(*this)( i,j ) = *ptr;                           // to access matrix coefficients,
}

//! @brief ... copy constructor
//! @param  rOther ... copied element
template<class T, int rows, int cols>
FullMatrix<T,rows,cols>::FullMatrix ( const FullMatrix<T,rows,cols>& rOther ): Eigen::Matrix<T,rows,cols>(rOther)
{
}

//! @brief ... constructor from Sparse matrices
//! @param  rOther ... copied element
template<class T, int rows, int cols>
FullMatrix<T,rows,cols>::FullMatrix ( const SparseMatrix<T>& rOther )
{
    size_t expectedSize =  sizeof(T) * rOther.GetNumColumns() * rOther.GetNumRows();
    size_t GB = 1e+9; // byte

    if (expectedSize > 5 * GB)
        throw MathException(__PRETTY_FUNCTION__, "Sparse matrix is too large for FullMatrix conversion. Approx. size of the FullMatrix would be " + std::to_string(expectedSize/GB) + "GB.");
    rOther.WriteEntriesToMatrix(*this);
}

//! @brief ... assignment constructor
//! @param  rOther ... copied element
template<class T, int rows, int cols>
template<typename OtherDerived>
FullMatrix<T,rows,cols>& FullMatrix<T,rows,cols>::operator=( const Eigen::MatrixBase <OtherDerived>& other)
{
	this->Eigen::Matrix<T,rows,cols>::operator=(other);
	return *this;
}

//! @brief ... assignment constructor for diagonal matrices
//! @param  rOther ... copied element
template<class T, int rows, int cols>
template<typename OtherDerived>
FullMatrix<T,rows,cols>& FullMatrix<T,rows,cols>::operator=( const Eigen::DiagonalBase< OtherDerived > &other)
{
	this->Eigen::Matrix<T,rows,cols>::operator=(other);
	return *this;
}

//! @brief ... assignment constructor for diagonal matrices
//! @param  rOther ... copied element
template<class T, int rows, int cols>
template<typename OtherDerived>
FullMatrix<T,rows,cols>& FullMatrix<T,rows,cols>::operator+=( const Eigen::DiagonalBase< OtherDerived > &other)
{
	this->Eigen::Matrix<T,rows,cols>::operator+=(other);
	return *this;
}

//! @brief ... resize matrix (everything is deleted, attention, the eigenroutine is not setting everything to zero)
//! @param rows ... number of rows
//! @param cols ... number of columns
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Resize ( int rRows, int rCols )
{
	if ( rRows==0 || rCols==0 )
		this->resize ( rRows,rCols );
	else
		this->setZero ( rRows, rCols );
}

#ifndef SWIG
//! @brief ... resize matrix (nothing is deleted)
//! @param rows ... number of rows
//! @param cols ... number of columns
//! routine is not exposed to python because of the different interface for vectors and matrices
//! if required, add the function manually in FullMatrix.i
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::ConservativeResize ( int rRows, int rCols )
{
	this->conservativeResize ( rRows,rCols );
}

//! @brief ... resize matrix (nothing is deleted)
//! @param rows ... number of rows
//! routine is not exposed to python because of the different interface for vectors and matrices
//! if required, add the function manually in FullMatrix.i
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::ConservativeResizeRows (int rRows )
{
	this->conservativeResize(rRows, Eigen::NoChange_t());
}

//! @brief ... resize matrix (nothing is deleted)
//! @param rows ... number of columns
//! routine is not exposed to python because of the different interface for vectors and matrices
//! if required, add the function manually in FullMatrix.i
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::ConservativeResizeCols (int rCols )
{
	this->conservativeResize(Eigen::NoChange_t(), rCols);
}
#endif

//! @brief ... add in place another matrix to this matrix and return a reference to this matrix
//! @param other ... other matrix
//! @return reference to this matrix
template<class T, int rows, int cols>
template<typename OtherDerived>
FullMatrix<T,rows,cols>& FullMatrix<T,rows,cols>::operator+=( const Eigen::MatrixBase <OtherDerived>& other)
{
	if ( this->rows() !=other.rows() || this->cols() !=other.cols() )
		throw MathException ( std::string ( "[FullMatrix::operator+=] Row or column number must be identical." ) );
	this->Eigen::Matrix<T,rows,cols>::operator+=(other);
	return *this;
}

//! @brief ... subtract in place another matrix from this matrix and return reference to this matrix
//! @param other ... other matrix
//! @return reference to this matrix
template<class T, int rows, int cols>
template<typename OtherDerived>
FullMatrix<T,rows,cols>& FullMatrix<T,rows,cols>::operator-=( const Eigen::MatrixBase <OtherDerived>& other)
{
	if ( this->rows() !=other.rows() || this->cols() !=other.cols() )
		throw MathException ( std::string ( "[FullMatrix::operator+=] Row or column number must be identical." ) );
	this->Eigen::Matrix<T,rows,cols>::operator-=(other);
	return *this;
}
//! @brief ... scale this matrix in place by a scalar factor and return a reference to this matrix
//! @param other ... scalar factor
//! @return reference to this matrix
template<class T, int rows, int cols>
FullMatrix<T,rows,cols>& FullMatrix<T,rows,cols>::operator*= ( const T &other)
{
	this->Eigen::Matrix<T,rows,cols>::operator*=(other);
	return *this;
}

//! @brief ... add a value to all entries of the matrix and return a reference to this matrix
//! @param other ... value to be added
//! @return reference to this matrix
template<class T, int rows, int cols>
FullMatrix<T,rows,cols>& FullMatrix<T,rows,cols>::operator+= ( const T &other)
{
	this->array()+=other;
	return *this;
}


//! @brief ... calculates the scalar product
//! @param other ... other matrix
//! @return reference to this matrix
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::Dot( const FullMatrix<T,rows,cols> &other )
{
	if ( this->rows() !=other.rows() || this->cols() !=1 || other.cols()!=1 )
		throw MathException ( std::string ( "[FullMatrix::Dot] Dot product requires number of rows to be identical (and only a single column)." ) );
	 return (this->transpose()*other).eval()(0,0);
}

//! @brief ... set the value of the (i,j)-th matrix entry
//! @param i ... row
//! @param j ... column
//! @param value ... value of the matrix entry
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::SetValue (int i, int j, const T& value )
{
	if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
	{
		throw MathException ( std::string ( "[FullMatrix::SetValue] Row or column number out of range." ) );
	}
	(*this)( i,j ) =value;
}

//! @brief ... set the value of the (i,j)-th matrix entry
//! @param i ... row
//! @param j ... column
//! @param value ... value of the matrix entry
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::SetValue (int i, const T& value )
{
	if ( i>=GetNumRows() || i<0 || GetNumColumns()!=1)
	{
		throw MathException ( std::string ( "[FullMatrix::SetValue] Row or column number out of range." ) );
	}
	(*this)( i,0 ) =value;
}

//! @brief ... get the value of the (i,j)-th vector entry
//! @param i ... row
//! @param j ... column
//! @return the value of the (i,j)-th matrix entry
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::GetValue (int i, int j ) const
{
	if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
	{
		throw MathException ( std::string ( "[FullMatrix::GetValue] Row  number out of range or not a vector." ) );
	}
	return (*this) ( i,j );
}

//! @brief ... get the value of the (i)-th vector entry
//! @param i ... row
//! @param j ... column
//! @return the value of the (i,j)-th matrix entry
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::GetValue (int i) const
{
	if ( i>=GetNumRows() || i<0 || GetNumColumns()!=1)
	{
		throw MathException ( std::string ( "[FullMatrix::GetValue] Row number out of range or not a vector." ) );
	}
	return (*this) ( i,0 );
}

//! @brief ... add a value to the (i,j)-th matrix entry
//! @param i ... row
//! @param j ... column
//! @param value ... value which is added to the (i,j)-th matrix entry
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::AddValue (int i, int j, const T& value )
{
	if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
		throw MathException ( std::string ( "[FullMatrix::AddValue] Row or column number out of range." ) );
	(*this) ( i,j ) +=value;
}

//! @brief ... subtract a value from the (i,j)-th matrix entry
//! @param i ... row
//! @param j ... column
//! @param value ... value which is subtracted from the (i,j)-th matrix entry
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::SubValue (int i, int j, const T& value )
{
	if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
		throw MathException ( std::string ( "[FullMatrix::SetValue] Row or column number out of range." ) );
	(*this) ( i,j )-=value;
}

//! @brief ... get number of rows
//! @return number of rows
template<class T, int rows, int cols>
inline int FullMatrix<T,rows,cols>::GetNumRows() const
{
	return this->rows();
}

//! @brief ... get number of columns
//! @return number of columns
template<class T, int rows, int cols>
inline int FullMatrix<T,rows,cols>::GetNumColumns() const
{
	return this->cols();
}


#ifdef ENABLE_SERIALIZATION
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Save ( const std::string &filename, std::string rType)const
{
    try
    {
	//transform to uppercase
	std::transform(rType.begin(), rType.end(), rType.begin(), toupper);

	// open file
	std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
	if (!ofs.is_open())
	{
	    throw MathException ( std::string ( "[NuTo::FullMatrix::Save] error opening file." ) );
	}

	// write data
	std::string typeIdString ( GetTypeId() );
	if (rType=="BINARY")
	{
	    boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
	    oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="XML")
	{
	    boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
	    oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="TEXT")
	{
	    boost::archive::text_oarchive ota ( ofs, std::ios::binary );
	    ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else
	{
	    throw MathException ( "[NuTo::FullMatrix::Save]File type not implemented." );
	}

	// close file
	ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
	std::string s ( std::string ( "[NuTo::FullMatrix::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
	throw MathException ( "[NuTo::FullMatrix::Save] Unhandled exception." );
    }
}

template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Restore ( const std::string &filename,  std::string rType)
{
    try
    {
	//transform to uppercase
	std::transform(rType.begin(), rType.end(), rType.begin(), toupper);

	// open file
	std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
	if (!ifs.is_open())
	{
	    throw MathException ( "[NuTo::FullMatrix::Restore] error opening file");
	}

	// read date
	std::string typeIdString;
	if (rType=="BINARY")
	{
	    boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
	    oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    if ( typeIdString!=GetTypeId() )
	    {
		throw MathException ( "[NuTo::FullMatrix::Restore]Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
	    }
	    oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="XML")
	{
	    boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
	    oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    if ( typeIdString!=GetTypeId() )
	    {
		throw MathException ( "[NuTo::FullMatrix::Restore]Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
	    }
	    oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="TEXT")
	{
	    boost::archive::text_iarchive ota ( ifs, std::ios::binary );
	    ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    if ( typeIdString!=GetTypeId() )
	    {
		throw MathException ( "[NuTo::FullMatrix::Restore]Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
	    }
	    ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else
	{
	    throw MathException ( "[NuTo::FullMatrix::Restore]File type not implemented" );
	}

	// close file
	ifs.close();
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
	throw MathException ( "[NuTo::FullMatrix::Restore] Unhandled exception." );
    }
}

#ifndef SWIG
    //! @brief serializes the class, this is the load routine
    //! @param ar         archive
    //! @param version    version
    template<class T, int rows, int cols>
    template<class Archive>
    void FullMatrix<T,rows,cols>::load(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    	std::cout << "start serialize FullMatrix (load)" << std::endl;
#endif
        std::vector<T> dataVec;
        int numRows, numColumns;
        ar & boost::serialization::make_nvp ("Matrix",boost::serialization::base_object< Matrix<T> > ( *this ) )
        & BOOST_SERIALIZATION_NVP(dataVec)
        & BOOST_SERIALIZATION_NVP(numRows)
        & BOOST_SERIALIZATION_NVP(numColumns);

        Resize(numRows,numColumns);
        memcpy ( this->data(),&(dataVec[0]),numRows * numColumns *sizeof ( T ) );
#ifdef DEBUG_SERIALIZATION
    	std::cout << "finish serialize FullMatrix (load)" << std::endl;
#endif
    }

    //! @brief serializes the class, this is the save routine
    //! @param ar         archive
    //! @param version    version
    template<class T, int rows, int cols>
    template<class Archive>
    void FullMatrix<T,rows,cols>::save(Archive & ar, const unsigned int version) const
    {
#ifdef DEBUG_SERIALIZATION
    	std::cout << "start serialize FullMatrix (save)" << std::endl;
#endif
        std::vector<T> dataVec ( GetNumRows() *GetNumColumns());
        memcpy ( & ( dataVec[0] ),this->data(),GetNumRows() *GetNumColumns() *sizeof ( T ) );
        int numRows = this->rows(),
                      numColumns = this->cols();
        ar & boost::serialization::make_nvp ("Matrix",boost::serialization::base_object< Matrix<T> > ( *this ) )
        & BOOST_SERIALIZATION_NVP(dataVec)
        & BOOST_SERIALIZATION_NVP(numRows)
        & BOOST_SERIALIZATION_NVP(numColumns);
#ifdef DEBUG_SERIALIZATION
    	std::cout << "finish serialize FullMatrix (save)" << std::endl;
#endif
    }

#endif
#endif //ENABLE_SERIALIZATION

//! @brief ... print info about the object
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Info() const
{
	Info ( 10,3 );
}

    //! @brief ... print info about the object
    //! @param width ... total width of each entry in the matrix when printed
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Info ( int width ) const
{
	Info ( width,3 );
}

//! @brief ... print info about the object
//! @param width ... total width of the each entry in the matrix when printed
//! @param precision ... precision (number of fractional digits)
//! @param rScientific ... switch for scientific notation
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Info ( int rWidth, int rPrecision, bool rScientific) const
{
	for ( int count=0; count<GetNumRows(); count++ )
	{
		for ( int count2=0; count2<GetNumColumns(); count2++ )
		{
			std::cout << this->Convert2String ( this->data() [count2*GetNumRows() +count], rScientific, rPrecision, rWidth) <<" ";
		}
		std::cout<<std::endl;
	}
}

//! @brief ... print info about the object into logger class
//! @param width ... total width of the each entry in the matrix when printed
//! @param precision ... precision (number of fractional digits)
//! @param rScientific ... switch for scientific notation
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Out(NuTo::Logger& rLogger, int rWidth, int rPrecision, bool rScientific)const
{
	for ( int count=0; count<GetNumRows(); count++ )
	{
		for ( int count2=0; count2<GetNumColumns(); count2++ )
		{
			rLogger << this->Convert2String ( this->data() [count2*GetNumRows() +count], rScientific, rPrecision, rWidth) <<" ";
		}
		rLogger << "\n";
	}
	rLogger << "\n";
}

//! @brief ... reads a matrix from a file
//! fileName_ ... file name
//! linesToSkip_ ... number of lines to skip at the beginning
//! delimiters_ ... delimiters between the entries in each line
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::ReadFromFile ( std::string fileName_, unsigned int linesToSkip_, std::string delimiters_ )
{
	try
	{
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

		// assign external file
		std::ifstream fileStream ( fileName_.c_str() );

		// check, if it could be opened
		if ( !fileStream.is_open() )
			throw MathException ( std::string ( "[FullMatrix::ReadFromFile]File could not be opened." ) );

		// get total number of lines
		fileStream.seekg (0, std::ios::beg);
		unsigned int totNumLines(0);
		std::string line;
		while (getline(fileStream, line))
		{
			totNumLines++;
			//std::cout << "line " << totNumLines << ": " << line << " (" << line.size() << ")" << std::endl;
		}

		if (linesToSkip_ > totNumLines)
		{
			throw MathException ("[FullMatrix::ReadFromFile] number of lines in file is smaller the the numbers of lines to be skipped.");
		}

		int numRows_ = totNumLines - linesToSkip_;
		//std::cout << "numRows_: " << numRows_ << std::endl;
		if (numRows_ > 0)
		{
			// reset stream
			fileStream.clear();
			fileStream.seekg (0, std::ios::beg);

			// read lines to be skipped
			for (unsigned int curLine_ =0; curLine_ < linesToSkip_; curLine_++)
			{
				getline(fileStream, line);
			}

			// read first line in order to obtain number of columns
			getline(fileStream, line);
			int numColumns_ ( 0 );
			boost::char_separator<char> sep ( delimiters_.c_str() );
			tokenizer tok ( line,sep );
			for ( tokenizer::iterator beg=tok.begin(); beg!=tok.end(); ++beg )
			{
				numColumns_++;
			}
			//std::cout << "numColumns_: " << numColumns_ << std::endl;

			// resize matrix
			this->resize ( numRows_,numColumns_ );

			// reset stream
			fileStream.clear();
			fileStream.seekg (0, std::ios::beg);

			// read lines to be skipped
			for (unsigned int curLine_ = 0; curLine_ < linesToSkip_; curLine_++)
			{
				getline(fileStream, line);
			}

			// read matrix
			for (int curRow_ = 0; curRow_ < numRows_; curRow_++)
			{
				getline (fileStream, line);
				tokenizer tok ( line,sep );

				int curColumn_ ( 0 );
				for ( tokenizer::iterator beg=tok.begin(); beg!=tok.end(); ++beg )
				{
					if ( curColumn_>numColumns_ )
						throw MathException ( std::string ( "[FullMatrix::ReadFromFile]Number of columns in line " ) +this->Int2String ( curRow_ + linesToSkip_ )
											  +std::string ( " which is " ) +this->Int2String ( curColumn_ ) +std::string ( " differs from the first line to be read(" )
											  +this->Int2String ( numColumns_ ) +std::string ( " columns)." ) );
					//std::cout << "curRow_: " << curRow_ << " curColumn_: " << curColumn_ << std::endl;
					(*this) ( curRow_,curColumn_ ) =Matrix<T>::ConvertFromString ( *beg );
					curColumn_++;
				}
			}
		}
		fileStream.close();
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
		throw MathException ( "[Matrix::ReadFromFile]Unhandled exception." );
	}
}

//! @brief ... reads a matrix from a file
//! @param  fileName_ ... file name
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::ReadFromFile ( std::string fileName_ )
{
	ReadFromFile ( fileName_,0," " );
}

//! @brief ... writes a matrix to a file
//! @param fileName_ ... file name
//! @param delimiter_ ... delimiters between the entries in each line
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::AppendToFile ( const std::string& fileName_, std::string delimiter_ ) const
{
	AppendToFile ( fileName_,delimiter_,"","" );
}

//! @brief ... append a matrix to a file
//! @param fileName_ ... file name
//! @param delimiter_ ... delimiters between the entries in each line
//! @param stringBefore ... String to write before the matrix (e.g. some comments)
//! @param stringAfter ... String to write after the matrix (e.g. some comments)
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::AppendToFile ( const std::string& fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter ) const
{
	try
	{
		// assign external file
		std::ofstream fileStream ( fileName_.c_str(),std::ios_base::out | std::ios_base::app );
		// check, if it could be opened
		if ( !fileStream.is_open() )
			throw MathException ( std::string ( "[FullMatrix::WriteToFile]File could not be opened." ) );
		fileStream << stringBefore;
		if ( stringBefore.length() >0 )
			fileStream << "\n";
		WriteToFile ( fileStream,delimiter_ );
		fileStream << stringAfter;
		fileStream.close();

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
		throw MathException ( "[Matrix::WriteToFile]Unhandled exception." );
	}
}

//! @brief ... writes a matrix to a file
//! @param fileName_ ... file name
//! @param delimiter_ ... delimiters between the entries in each line
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::WriteToFile ( const std::string& fileName_, std::string delimiter_ ) const
{
	WriteToFile ( fileName_,delimiter_,"","" );
}


//! @brief ... writes a matrix to a file
//! @param fileName_ ... file name
//! @param delimiter_ ... delimiters between the entries in each line
//! @param stringBefore ... String to write before the matrix (e.g. some comments)
//! @param stringAfter ... String to write after the matrix (e.g. some comments)
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::WriteToFile ( const std::string& fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter ) const
{
	try
	{
		// assign external file
		std::ofstream fileStream ( fileName_.c_str() );
		// check, if it could be opened
		if ( !fileStream.is_open() )
			throw MathException ( std::string ( "[FullMatrix::WriteToFile]File could not be opened." ) );
		fileStream << stringBefore;
		if ( stringBefore.length() >0 )
			fileStream << "\n";
		WriteToFile ( fileStream,delimiter_ );
		fileStream << stringAfter;
		fileStream.close();
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
		throw MathException ( "[Matrix::WriteToFile]Unhandled exception." );
	}
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
template<class T, int rows, int cols>
std::string FullMatrix<T,rows,cols>::GetTypeId()const
{
	return std::string("FullMatrix")+GetTypeIdBaseDataType<T>();
}


//! @brief ... converts any FullMatrix to a DoubleFullMatrix
//! @return    converted matrix
template<class T, int rows, int cols>
FullMatrix<double,rows,cols> FullMatrix<T,rows,cols>::Convert2double()
{
	FullMatrix<double,rows,cols> doubleMatrix(GetNumRows(),GetNumColumns());
	for (int count=0; count<GetNumColumns(); count++)
		for (int count2=0; count2<GetNumRows(); count2++)
			doubleMatrix(count2,count) = (double)(*this)(count2,count);

	return doubleMatrix;
}

//! @brief ... converts any FullMatrix to a IntFullMatrix
//! @return    converted matrix
template<class T, int rows, int cols>
FullMatrix<int,rows,cols> FullMatrix<T,rows,cols>::Convert2int()
{
	FullMatrix<int,rows,cols> intMatrix(GetNumRows(),GetNumColumns());
	for (int count=0; count<GetNumColumns(); count++)
		for (int count2=0; count2<GetNumRows(); count2++)
			intMatrix(count2,count) = (int)(*this)(count2,count);

	return intMatrix;
}

//! @brief solves a system of linear equations with a symmetric positive definite coefficient matrix using Cholesky factorization (LAPACK)
//! @param rRHS ... matrix of right-hand-side vectors
//! @param rSolution ... solution matrix
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::SolveCholeskyLapack(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRHS,
		FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSolution) const
{
	return this->SolveCholeskyLapack(rRHS, rSolution, params<T,rows,cols>());
}


//! @brief this is the norm member function that is called for anything else than double
template<typename T, int rows, int cols>
template<typename T1, int rows1, int cols1>
void FullMatrix<T,rows,cols>::SolveCholeskyLapack(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRHS,
		FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSolution, params<T1, rows1, cols1>)const
{
	throw MathException("[FullMatrix::SolveCholeskyLapack] only implemented for double values.");
}

//! @brief this is the norm member function that is called for doubles
template<typename T, int rows, int cols>
template<int rows1, int cols1>
void FullMatrix<T,rows,cols>::SolveCholeskyLapack(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRHS,
		FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSolution, params<double, rows1, cols1>)const
{
#ifdef ENABLE_MKL
	// check matrix
	if(this->GetNumColumns() != this->GetNumRows())
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid shape of the coefficient matrix.");
	}
	if(this->GetNumColumns() < 1)
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid dimension of the coefficient matrix.");
	}
	// check right-hand-side vectors
	if(this->GetNumRows() != rRHS.GetNumRows())
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid number of rows in right-hand-side matrix.");
	}
	if(rRHS.GetNumColumns() < 1)
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid number of columns in right-hand-side matrix.");
	}

	// copy coefficient matrix and perform cholesky factorization
	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> choleskyMatrix(*this);
	int dimCholeskyMatrix = choleskyMatrix.GetNumColumns();
	char uplo('L');
	int info(0);
	dpotrf_(&uplo, &dimCholeskyMatrix, choleskyMatrix.mEigenMatrix.data(), &dimCholeskyMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::SolveCholeskyLapack] Error calculating Cholesky factorization.");
		if(info < 1)
		{
			info *= -1;
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		}
		else
		{
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The leading minor of order "+infoStream.str()+"is not positive definite.";
		}
		throw MathException(errorMessage);
	}

	// copy right-hand-side matrix and solve linear system of equations using factorized choolesky matrix
	rSolution = rRHS;
	int numRHSVectors = rSolution.GetNumColumns();
	dpotrs_(&uplo, &dimCholeskyMatrix, &numRHSVectors, choleskyMatrix.mEigenMatrix.data(), &dimCholeskyMatrix, rSolution.mEigenMatrix.data(), &dimCholeskyMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::SolveCholeskyLapack] Error solving linear system of equations.");
		std::stringstream infoStream;
		infoStream << -1*info;
		errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		throw MathException(errorMessage);
	}
#else //ENABLE_MKL
	throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] lapack package not enabled.");
#endif //ENABLE_MKL
}

//! @brief calculates the inverse of a symmetric positive definite matrix using Cholesky factorization (LAPACK)
//! @param rInverse ... inverse matrix
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::InverseCholeskyLapack(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInverse) const
{
	return this->InverseCholeskyLapack(rInverse, params<T,rows,cols>());
}

//! @brief this is the InverseCholeskyLapack member function that is called for anything else than double
template<typename T, int rows, int cols>
template<typename T1, int rows1, int cols1>
void FullMatrix<T,rows,cols>::InverseCholeskyLapack(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInverse, params<T1, rows1, cols1>)const
{
	throw MathException("[FullMatrix::InverseCholeskyLapack] only implemented for double values.");
}

//! @brief this is the InverseCholeskyLapack member function that is called for doubles
template<typename T, int rows, int cols>
template<int rows1, int cols1>
void FullMatrix<T,rows,cols>::InverseCholeskyLapack(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInverse, params<double, rows1, cols1>)const
{
#ifdef ENABLE_MKL
	// check matrix
	if(this->GetNumColumns() != this->GetNumRows())
	{
		throw MathException("[NuTo::FullMatrix::InverseCholeskyLapack] invalid shape of the coefficient matrix.");
	}
	if(this->GetNumColumns() < 1)
	{
		throw MathException("[NuTo::FullMatrix::InverseCholeskyLapack] invalid dimension of the coefficient matrix.");
	}

	// copy matrix and perform cholesky factorization
	rInverse = *this;
	int dimMatrix = rInverse.GetNumColumns();
	char uplo('L');
	int info(0);
	dpotrf_(&uplo, &dimMatrix, rInverse.data(), &dimMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::InverseCholeskyLapack] Error calculating Cholesky factorization.");
		if(info < 1)
		{
			info *= -1;
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		}
		else
		{
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The leading minor of order "+infoStream.str()+"is not positive definite.";
		}
		throw MathException(errorMessage);
	}

	// calculate inverse matrix
	dpotri_(&uplo, &dimMatrix, rInverse.mEigenMatrix.data(), &dimMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::InverseCholeskyLapack] Error calculating inverse matrix.");
		if(info < 1)
		{
			info *= -1;
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		}
		else
		{
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The ("+infoStream.str()+","+infoStream.str()+")-th element of the factor U or L is zero, and the inverse could not be computed.";
		}
		throw MathException(errorMessage);
	}

	// copy values from lower triangle to upper triangle
	double* data = rInverse.data();
	for(int col = 1; col < dimMatrix; col++)
	{
		for(int row = 0; row < col; row++)
		{
			data[col * dimMatrix + row] = data[row * dimMatrix + col];
		}
	}
#else //ENABLE_MKL
	throw MathException("[NuTo::FullMatrix::InverseCholeskyLapack] lapack package not enabled.");
#endif //ENABLE_MKL
}


//! @brief ... calculates the norm of this matrix, i.e. for vectors the Euclidean norm
//! @return norm of this matrix
template<class T, int rows, int cols>
FullMatrix<T,rows,cols>  FullMatrix<T,rows,cols>::Inverse() const
{
	return this->Inverse(params<T,rows,cols>());
}

//! @brief this is the norm member function that is called for anything else than double
template<typename T, int rows, int cols>
template<typename T1, int rows1, int cols1>
FullMatrix<T,rows,cols> FullMatrix<T,rows,cols>::Inverse(params<T1, rows1, cols1>)const
{
	throw MathException("[FullMatrix::Inverse] only implemented for double values.");
}

//! @brief this is the norm member function that is called for doubles
template<typename T, int rows, int cols>
template<int rows1, int cols1>
FullMatrix<T,rows,cols> FullMatrix<T,rows,cols>::Inverse(params<double, rows1, cols1>)const
{
	return this->inverse().eval();
}

//! @brief calculates the eigenvalues
//! @param rEigenValues ... eigenvalues
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::EigenValuesSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues) const
{
	return this->EigenValuesSymmetric(rEigenValues, params<T,rows,cols>());
}

//! @brief this is the EigenValuesSymmetric member function that is called for anything else than double
template<typename T, int rows, int cols>
template<typename T1, int rows1, int cols1>
void FullMatrix<T,rows,cols>::EigenValuesSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, params<T1, rows1, cols1>)const
{
	throw MathException("[FullMatrix::EigenValuesSymmetric] only implemented for double values.");
}

//! @brief this is the EigenValuesSymmetric member function that is called for doubles
template<typename T, int rows, int cols>
template<int rows1, int cols1>
void FullMatrix<T,rows,cols>::EigenValuesSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, params<double, rows1, cols1>)const
{
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> mySolver((*this),false);
	rEigenValues = FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>(mySolver.eigenvalues());
}


//! @brief calculates the eigenvectors
//! @param rEigenVectors ... eigenvectors
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::EigenVectorsSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors)const
{
	return this->EigenVectorsSymmetric(rEigenValues, rEigenVectors, params<T,rows,cols>());
}


//! @brief this is the EigenVectorsSymmetric member function that is called for anything else than double
template<typename T, int rows, int cols>
template<typename T1, int rows1, int cols1>
void FullMatrix<T,rows,cols>::EigenVectorsSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<T1, rows1, cols1>)const
{
	throw MathException("[FullMatrix::EigenVectorsSymmetric] only implemented for double values.");
}

//! @brief this is the EigenVectorsSymmetric member function that is called for doubles
template<typename T, int rows, int cols>
template<int rows1, int cols1>
void FullMatrix<T,rows,cols>::EigenVectorsSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<double, rows1, cols1>)const
{
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> mySolver((*this));
	rEigenValues  = mySolver.eigenvalues();
	rEigenVectors = mySolver.eigenvectors();
}

//! @brief calculates the eigenvectors
//! @param rEigenVectors ... eigenvectors
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::GeneralizedEigenVectorsSymmetric(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rM, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors)const
{
	return this->GeneralizedEigenVectorsSymmetric(rM, rEigenValues, rEigenVectors, params<T,rows,cols>());
}


//! @brief this is the EigenVectorsSymmetric member function that is called for anything else than double
template<typename T, int rows, int cols>
template<typename T1, int rows1, int cols1>
void FullMatrix<T,rows,cols>::GeneralizedEigenVectorsSymmetric(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rM,FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<T1, rows1, cols1>)const
{
	throw MathException("[FullMatrix::GeneralizedEigenVectorsSymmetric] only implemented for double values.");
}

//! @brief this is the EigenVectorsSymmetric member function that is called for doubles
template<typename T, int rows, int cols>
template<int rows1, int cols1>
void FullMatrix<T,rows,cols>::GeneralizedEigenVectorsSymmetric(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rM, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<double, rows1, cols1>)const
{
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> mySolver((*this),rM);
	rEigenValues  = mySolver.eigenvalues();
	rEigenVectors = mySolver.eigenvectors();
}

//! @brief ... imports a matrix from a SLang ASCII file
//! @param fileName ... file name
//template<class T, int rows, int cols>
//void FullMatrix<T,rows,cols>::ImportFromSLangText ( const char* fileName )
//{
//	throw MathException("to be done");
//}

//! @brief ... imports Matrix from a Vtk ASCII File
//! @param fileName ... file name
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::ImportFromVtkASCIIFile(const char* rFileName)
{
    using namespace boost::spirit::classic;

    // open file
    std::ifstream file(rFileName, std::ios::in);
    if (file.is_open() == false)
    {
        throw MathException("[FullMatrix::ImportFromVtkASCIIFile] error opening file.");
    }
    //skip read header, is done in StructureGrid.cpp
    // read first four lines
    int numEntries(0);

    std::string line;
    for (int count=0;count<7;count++)
    {
        getline (file, line);
    }

    // read number of entries
    getline (file, line);
    if (parse(line.c_str(),("POINT_DATA ">> int_p[assign_a(numEntries)] >>  *space_p)).full == false)
           {
               throw MathException("[Matrix::importFromVtkASCIIFile]error reading number of entries.");
           }
    // read data type
    getline (file, line);
    // read empty line
    getline (file, line);

  // read entries
    std::vector<T> imageValues;
    T value;
    imageValues.reserve(numEntries);

    while(getline(file,line))
    {
    	std::istringstream iss(line);
    	while(iss >> value)
    	{
			imageValues.push_back(value);
    	}
    }
    for (int count = 0; count<numEntries; count++)
    {
        (*this)(count, 0) = imageValues[count];
    }
    // close file
   file.close();
}

//! @brief performs a dyadic operator on all matrix entries and the entries of another matrix
//! @param rDOperator ... Dyadic Operator
//! @param otherMatrix ... other matrix
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Map ( NuTo::DyadicOperator<T>* rDOperator, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &otherMatrix )
{
	throw MathException ( "FullMatrix::Map - To be implemented." );
}

//! @brief performs a monadic operator on all matrix entries
//! @param rMOperator ... Monadic Operator
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Map ( const NuTo::MonadicOperator<T>* rMOperator )
{
	for ( int count=0; count<GetNumColumns() *GetNumRows(); count++ )
		this->data() [count] = rMOperator->Evaluate ( this->data() [count] );
}

//! @brief performs a dyadic operator on all matrix entries with another given value
//! @param rDOperator ... Dyadic Operator
//! @param rValue ... value
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::Map ( const NuTo::DyadicOperator<T>* rDOperator, const T& rValue )
{
	for ( int count=0; count<GetNumColumns() *GetNumRows(); count++ )
		this->data() [count] = rDOperator->Evaluate ( this->data() [count],rValue );
}

//! @brief returns the maximum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
//! @param rRowOutput ... row
//! @param rColumnOutput ... column
//! @return maximum value of the matrix
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::Max ( int& rRowOutput, int& rColumnOutput ) const
{
	if ( GetNumColumns() ==0 || GetNumRows() ==0 )
		throw MathException ( "FullMatrix::Max - Maximum for matrix with zero entries cannot be calculated." );

	return this->maxCoeff ( &rRowOutput,&rColumnOutput );
}

//! @brief returns the maximum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
//! @return maximum value of the matrix
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::Max() const
{
	if ( GetNumColumns() ==0 || GetNumRows() ==0 )
		throw MathException ( "FullMatrix::Max - Maximum for matrix with zero entries cannot be calculated." );
	return this->maxCoeff();
}

//! @brief returns the minimum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
//! @param rRowOutput ... row
//! @param rColumnOutput ...  column
//! @return minimum value of the matrix
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::Min ( int& rRowOutput, int& rColumnOutput ) const
{
	if ( GetNumColumns() ==0 || GetNumRows() ==0 )
		throw MathException ( "FullMatrix::Min - Minimum for matrix with zero entries cannot be calculated." );

	return this->minCoeff ( &rRowOutput,&rColumnOutput );
}

//! @brief returns the minimum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
//! @return minimum value of the matrix
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::Min() const
{
	if ( GetNumColumns() ==0 || GetNumRows() ==0 )
		throw MathException ( "FullMatrix::Min - Minimum for matrix with zero entries cannot be calculated." );

	return this->minCoeff();
}

//! @brief elementwise absolute value of the matrix
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::Abs() const
{
	return this->array().abs().matrix();
}

//! @brief calcualte the transpose of this matrix
//! @return transpose of the matrix
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::Trans() const
{
	return this->transpose();
}

//! @brief ... extract a block (submatrix) from this matrix
//! @param rI ... start row
//! @param rJ ... start column
//! @param rRows ... number of rows of the submatrix
//! @param rCols ... number of columns of the submatrix
//! @return a submatrix with size (rRows,rCols) extracted at (rI,rJ)-entry of this matrix
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>  FullMatrix<T,rows,cols>::GetBlock ( int rI, int rJ, int rRows, int rCols ) const
{
	if ( rI+rRows>this->rows() )
		throw MathException ( "FullMatrix::GetBlock - rows out of Dimension." );
	if ( rJ+rCols>this->cols() )
		throw MathException ( "FullMatrix::GetBlock - columns out of Dimension." );
	if ( rI<0 )
		throw MathException ( "FullMatrix::GetBlock - row should not be negative." );
	if ( rJ<0 )
		throw MathException ( "FullMatrix::GetBlock - column should not be negative." );
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->block ( rI,rJ,rRows,rCols ) );
}

//! @brief ... set a block (submatrix) in this Matrix
//! @param rI ... start row
//! @param rJ ... start column
//! @param rBlock ... submatrix
template<class T, int rows, int cols>
void  FullMatrix<T,rows,cols>::SetBlock ( int rI, int rJ, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock )
{
	if ( rI+rBlock.rows() >this->rows() )
		throw MathException ( "FullMatrix::SetBlock - rows out of Dimension." );
	if ( rJ+rBlock.cols() >this->cols() )
		throw MathException ( "FullMatrix::SetBlock - columns out of Dimension." );
	if ( rI<0 )
		throw MathException ( "FullMatrix::SetBlock - row should not be negative." );
	if ( rJ<0 )
		throw MathException ( "FullMatrix::SetBlock - column should not be negative." );
	this->block ( rI,rJ,rBlock.rows(),rBlock.cols() ) = rBlock;
}

//! @brief ... add a block (submatrix) in this Matrix
//! @param rI ... start row
//! @param rJ ... start column
//! @param rBlock ... submatrix
template<class T, int rows, int cols>
void  FullMatrix<T,rows,cols>::AddBlock ( int rI, int rJ, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock )
{
	if ( rI+rBlock.rows() >this->rows() )
		throw MathException ( "FullMatrix::SetBlock - rows out of Dimension." );
	if ( rJ+rBlock.cols() >this->cols() )
		throw MathException ( "FullMatrix::SetBlock - columns out of Dimension." );
	if ( rI<0 )
		throw MathException ( "FullMatrix::SetBlock - row should not be negative." );
	if ( rJ<0 )
		throw MathException ( "FullMatrix::SetBlock - column should not be negative." );
	this->block ( rI,rJ,rBlock.rows(),rBlock.cols() ) += rBlock;
}

//! @brief ... extract a row from the matrix
//! @param rI ... row
//! @return row as full matrix
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>  FullMatrix<T,rows,cols>::GetRow ( int rI ) const
{
	if ( rI>=this->rows() )
		throw MathException ( "FullMatrix::GetRow - row out of Dimension." );
	if ( rI<0 )
		throw MathException ( "FullMatrix::GetRow - row should not be negative." );
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->row ( rI ) );
}

//! @brief ... set a row in this FullMatrix
//! @param rI ... row
//! @param rBlock ... new row intries
template<class T, int rows, int cols>
void  FullMatrix<T,rows,cols>::SetRow ( int rI, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock )
{
	if ( rBlock.rows() !=1 )
		throw MathException ( "FullMatrix::SetRow - Expect a Matrix with a single row as Input." );
	if ( rI<0 )
		throw MathException ( "FullMatrix::SetRow - row should not be negative." );
	if ( rBlock.cols() !=this->cols() )
		throw MathException ( "FullMatrix::SetRow - number of columns for both matrices must be identical." );
	this->row ( rI ) = rBlock;
}

//! @brief ... extract a column from this matrix
//! @param rI ... column
//! @return column entries as full matrix
template<class T, int rows, int cols>
FullVector<T, Eigen::Dynamic>  FullMatrix<T,rows,cols>::GetColumn ( int rI ) const
{
	if ( rI>=this->cols() )
		throw MathException ( "FullMatrix::GetColumn - column out of Dimension." );
	if ( rI<0 )
		throw MathException ( "FullMatrix::GetColumn - column should not be negative." );
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->col ( rI ) );
}

//! @brief ... set a column in this matrix
//! @param rI ... column
//! @param rBlock ... new column entries
template<class T, int rows, int cols>
void  FullMatrix<T,rows,cols>::SetColumn ( int rI, const FullVector<T, Eigen::Dynamic>& rBlock )
{
	if ( rBlock.cols() !=1 )
		throw MathException ( "FullMatrix::SetColumn - Expect a Matrix with a single column as Input." );
	if ( rI<0 )
		throw MathException ( "FullMatrix::SetColumn - column should not be negative." );
	if ( rBlock.rows() !=this->rows() )
		throw MathException ( "FullMatrix::SetColumn - number of rows for both matrices must be identical." );
	this->col ( rI ) = rBlock;
}

//! @brief ... appends columns to this matrix
//! @param rBlock ... matrix storing the columns
template<class T, int rows, int cols>
void  FullMatrix<T,rows,cols>::AppendColumns (const FullVector<T, Eigen::Dynamic>& rBlock )
{
	if ( rBlock.rows() !=this->rows() )
		throw MathException ( "FullMatrix::AppendColumns - number of rows for both matrices must be identical." );
	Eigen::Matrix<T,rows,cols> tmpMatrix ( this->rows(),this->cols() +rBlock.cols() );
	tmpMatrix.block ( 0,0,this->rows(),this->cols() ) = *this;
	tmpMatrix.block ( 0,this->cols(),rBlock.rows(),rBlock.cols() ) = rBlock;
	*this = tmpMatrix;
}

//! @brief ... appends columns to this matrix
//! @param rBlock ... matrix storing the columns
template<class T, int rows, int cols>
void  FullMatrix<T,rows,cols>::AppendColumns (const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock )
{
	if ( rBlock.rows() !=this->rows() )
		throw MathException ( "FullMatrix::AppendColumns - number of rows for both matrices must be identical." );
	Eigen::Matrix<T,rows,cols> tmpMatrix ( this->rows(),this->cols() +rBlock.cols() );
	tmpMatrix.block ( 0,0,this->rows(),this->cols() ) = *this;
	tmpMatrix.block ( 0,this->cols(),rBlock.rows(),rBlock.cols() ) = rBlock;
	*this = tmpMatrix;
}

//! @brief ... appends rows to this matrix
//! @param rBlock ... matrix storing the rows
template<class T, int rows, int cols>
void  FullMatrix<T,rows,cols>::AppendRows (const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock )
{
	if ( rBlock.cols() !=this->cols() )
		throw MathException ( "FullMatrix::AppendRows - number of columns for both matrices must be identical." );
	Eigen::Matrix<T,rows,cols> tmpMatrix ( this->rows() +rBlock.rows(),this->cols() );
	tmpMatrix.block ( 0,0,this->rows(),this->cols() ) = (*this);
	tmpMatrix.block ( this->rows(),0,rBlock.rows(),rBlock.cols() ) = rBlock;
	*this = tmpMatrix;
}

//! @brief ... calculates the norm of this matrix, i.e. for vectors the Euclidean norm
//! @return norm of this matrix
template<class T, int rows, int cols>
double  FullMatrix<T,rows,cols>::Norm() const
{
	return this->Norm(params<T,rows,cols>());
}

//! @brief this is the norm member function that is called for anything else than double
template<typename T, int rows, int cols>
template<typename T1, int rows1, int cols1>
double FullMatrix<T,rows,cols>::Norm(params<T1, rows1, cols1>)const
{
	throw MathException("[FullMatrix::Norm] only implemented for real values.");
}

//! @brief this is the norm member function that is called for doubles
template<typename T, int rows, int cols>
template<int rows1, int cols1>
double FullMatrix<T,rows,cols>::Norm(params<double, rows1, cols1>)const
{
	return this->norm();
}

//! @brief ... computes the sum of all coefficients
//! @return sum of the entries for each column
template<class T, int rows, int cols>
T FullMatrix<T,rows,cols>::Sum() const
{
	return this->sum();
}

//! @brief ... sets small matrix entries to zero - e.g. for output formating
//! @param rTolerance ... cut-off tolerance
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::SetSmallEntriesZero(double rTolerance)
{
    for (int j = 0; j < this->cols(); ++j)
        for (int i = 0; i < this->rows(); ++i)
            if (std::abs((*this)(i, j)) < rTolerance)
                (*this)(i, j) = 0;
}

//! @brief ... calculates the sum of the entries for each column
//! @return sum of the entries for each column
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::ColumnwiseSum() const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->colwise().sum() );
}

//! @brief ... calculates the sum of the entries for each row
//! @return sum of the entries for each row
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::RowwiseSum() const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->rowwise().sum() );
}

//! @brief ... calculates the minimum of the entries for each column
//! @return the minimum of the entries for each column
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::ColumnwiseMinCoeff() const
{
	return this->colwise().minCoeff();
}

//! @brief ... calculates the minimum of the entries for each row
//! @return the minimum of the entries for each row
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::RowwiseMinCoeff() const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->rowwise().minCoeff() );
}

//! @brief ... calculates the maximum of the entries for each column
//! @return the maximum of the entries for each column
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::ColumnwiseMaxCoeff() const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->colwise().maxCoeff() );
}

//! @brief ... calculates the maximum of the entries for each row
//! @return the maximum of the entries for each row
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::RowwiseMaxCoeff() const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->rowwise().maxCoeff() );
}

//! @brief ... coefficient wise multiplication of this matrix with another matrix
//! @param other ... other matrix
//! @return a matrix which is obtained by a coefficient wise multiplication of this matrix with another matrix
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::ElementwiseMul ( const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &other ) const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( (this->array() *other.array()).matrix() );
}

//! @brief ... coefficient wise division of this matrix by another matrix
//! @param other ... other matrix
//! @return a matrix which is obtained by a coefficient wise division of this matrix by another matrix
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::ElementwiseDiv ( const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &other ) const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( (this->array() /other.array()).matrix() );
}

//! @brief ... coefficient wise reciprocal
//! @return a matrix which is obtained by coefficient wise reciprocal of this matrix
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::ElementwiseInverse() const
{
	return FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ( this->array().inverse().matrix());
}

//! @brief ... sorts the rows of a matrix based on the entries in column rCol
template<class T, int rows, int cols>
FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> FullMatrix<T,rows,cols>::SortRow(int rCol)
{
	if (rCol<0 || rCol>=this->GetNumColumns())
		throw MathException("[FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>::SortRow] incorrect row");

	struct mySortClassRow {
		const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> *mat;
		int colnum;
		mySortClassRow(const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &m, int keycol)  : mat(&m), colnum(keycol){};
		bool operator()(const int &lhs, const int &rhs) const
		{
			return (*mat)(lhs,colnum) < (*mat)(rhs,colnum);
		}
	};

	std::vector<unsigned int> theRows(this->GetNumRows());
	for (unsigned int count=0; count<theRows.size(); count++)
	{
		theRows[count]=count;
	}

	//sort the index vector according to the entries in the matrix
	std::sort(theRows.begin(),theRows.end(),mySortClassRow(*this,rCol));

	//rearrange the full matrix rows
	//for simplicity, just copy the matrix
	NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> myOrigData(*this);
	for (unsigned int count=0; count<theRows.size(); count++)
	{
		for ( int count2=0; count2<this->GetNumColumns(); count2++)
		{
			(*this)(count,count2) = myOrigData(theRows[count],count2);
		}
	}
	return *this;
}

//! @brief ... writes a matrix to a file
//! @param fileStream ... output file stream
//! @param delimiter_ ... delimiters between the entries in each line
template<class T, int rows, int cols>
void FullMatrix<T,rows,cols>::WriteToFile ( std::ofstream& fileStream, std::string& delimiter_ ) const
{
	// go through the matrix and store the values
	fileStream.setf ( std::ios::scientific, std::ios::floatfield );
	for ( int count=0; count<GetNumRows(); count++ )
	{
		for ( int count2=0; count2<GetNumColumns(); count2++ )
		{
			fileStream    <<     this->Convert2String ( this->data() [count+count2*GetNumRows() ],true,12,18 );
			if ( count2!=GetNumColumns()-1 )
				fileStream    << delimiter_;
		}
		if ( count<GetNumRows()-1 )
			fileStream    <<  "\n";
	}
}


} //NAMESPACE NUTO

#endif // FULL_MATRIX_H
