// $Id$

#ifndef FULL_MATRIX_H
#define FULL_MATRIX_H

#include <iomanip>  //for setw
#include <fstream>  //for file acces
#include <vector>
#include <string>
#include <iostream>

#include <boost/tokenizer.hpp>

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
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

#include "nuto/math/Matrix.h"
#include "nuto/math/NuToMath.h"
#include "nuto/math/MathException.h"
#include "nuto/math/fortran_routines.h"
#include "nuto/math/SparseMatrix.h"
#include <eigen2/Eigen/Core>

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date July 2009
//! @brief ... class for full matrices derived from the abstract base class Matrix
template <class T>
class FullMatrix : public Matrix<T>

{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif

public:
    //! @brief ... default constructor
    FullMatrix<T>()
    {
        mEigenMatrix.resize ( 0,0 );
    }

    //! @brief ... constructor
    //! @param rNumRows ... number of columns
    //! @param rNumColumns ... number of rows
    FullMatrix<T> ( int rNumRows, int rNumColumns)
    {
        if ( rNumRows*rNumColumns>0 )
            mEigenMatrix.setZero ( rNumRows,rNumColumns );
        else
            mEigenMatrix.resize ( rNumRows,rNumColumns );
    }

#ifndef SWIG
    //! @brief ... constructor
    //! @param rEigenMatrix ... other matrix
    FullMatrix<T> ( const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &rEigenMatrix )
    {
        mEigenMatrix = rEigenMatrix;
    }
#endif

    //! @brief ... constructur
    //! @param numRows_ ... number of columns
    //! @param numColumns_ ... number of rows
    //! @param entries_ ... vector containing the matrix in column-major orientation
    FullMatrix<T> ( int numRows_, int numColumns_, const std::vector<T>& entries_ )
    {
        mEigenMatrix.resize ( numRows_,numColumns_ );
        const T *ptr = &entries_[0];
        for ( int j=0; j<mEigenMatrix.cols(); ++j )              // loop over columns
            for ( int i=0; i<mEigenMatrix.rows(); ++i, ptr++ )   // loop over rows
                mEigenMatrix ( i,j ) = *ptr;                     // to access matrix coefficients,
    }

    //! @brief ... constructor
    //! @brief ... creates a FullVector(means a FullMatrix size(m,1) )
    //! @param entries_ ... vector containing the matrix in column-major orientation
    FullMatrix<T> ( const std::vector<T>& entries_ )
    {
        mEigenMatrix.resize ( entries_.size(),1 );
        const T *ptr = &entries_[0];

		for ( int i=0; i<mEigenMatrix.rows(); i++, ptr++ )   	// loop over rows
			mEigenMatrix ( i,0 ) = *ptr;                     	// to access matrix coefficients,
    }

   //! @brief ... copy constructor
    //! @param  rOther ... copied element
    FullMatrix<T> ( const FullMatrix<T>& rOther )
    {
        mEigenMatrix = rOther.mEigenMatrix;
    }

    //! @brief ... resize matrix (everything is deleted)
    //! @param rows ... number of columns
    //! @param cols ... number of rows
    void Resize ( int rows, int cols )
    {
    	if ( rows*cols>0 )
            mEigenMatrix.setZero ( rows, cols );
        else
            mEigenMatrix.resize ( rows,cols );
    }

//    template <class F>
//    friend std::ostream& operator<<(std::ostream& os, const FullMatrix<F>& r);

    FullMatrix<T> operator+ ( const FullMatrix<T> &other ) const
    {
        if ( mEigenMatrix.rows() !=other.mEigenMatrix.rows() || mEigenMatrix.cols() !=other.mEigenMatrix.cols() )
            throw MathException ( std::string ( "[FullMatrix::operator+] Row or column number must be identical." ) );
        return FullMatrix<T> ( mEigenMatrix+other.mEigenMatrix );
    }

    FullMatrix<T> operator- ( const FullMatrix<T> &other ) const
    {
        if ( mEigenMatrix.rows() !=other.mEigenMatrix.rows() || mEigenMatrix.cols() !=other.mEigenMatrix.cols() )
            throw MathException ( std::string ( "[FullMatrix::operator-] Row or column number must be identical." ) );
        return FullMatrix<T> ( mEigenMatrix-other.mEigenMatrix );
    }

    FullMatrix<T> operator* ( const FullMatrix<T> &other ) const
    {
    	if ( mEigenMatrix.cols() !=other.mEigenMatrix.rows() )
        	throw MathException ( std::string ( "[FullMatrix::operator*] Number of columns of the first matrix must be identical to the number of rows of the second matrix." ) );

    	return FullMatrix<T> ( mEigenMatrix*other.mEigenMatrix );
    }

    FullMatrix<T> operator* ( const T &other ) const
    {
        return FullMatrix<T> ( mEigenMatrix*other );
    }

    FullMatrix<T>& operator+= ( const FullMatrix<T> &other )
    {
        if ( mEigenMatrix.rows() !=other.mEigenMatrix.rows() || mEigenMatrix.cols() !=other.mEigenMatrix.cols() )
            throw MathException ( std::string ( "[FullMatrix::operator+=] Row or column number must be identical." ) );
        mEigenMatrix+=other.mEigenMatrix;
        return *this;
    }

    FullMatrix<T>& operator-= ( const FullMatrix<T> &other )
    {
        if ( mEigenMatrix.rows() !=other.mEigenMatrix.rows() || mEigenMatrix.cols() !=other.mEigenMatrix.cols() )
            throw MathException ( std::string ( "[FullMatrix::operator-=] Row or column number must be identical." ) );
        mEigenMatrix-=other.mEigenMatrix;
        return *this;
    }
    
    FullMatrix<T> operator*= ( const T &other )
    {
        mEigenMatrix*=other;
        return *this;
    }


    T& operator() ( int i, int j )
    {
        if ( i>=mEigenMatrix.rows() || j>=mEigenMatrix.cols() || i<0 || j<0 )
            throw MathException ( std::string ( "[FullMatrix::operator()] Row or column number out of range." ) );
        return  mEigenMatrix ( i,j );
    }

#ifndef SWIG
    const T operator() ( int i, int j ) const
    {
        if ( i>=mEigenMatrix.rows() || j>=mEigenMatrix.cols() || i<0 || j<0 )
            throw MathException ( std::string ( "[FullMatrix::operator()] Row or column number out of range." ) );
        return  mEigenMatrix ( i,j );
    }
#endif

    void SetValue (int i, int j, const T& value )
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
        {
            throw MathException ( std::string ( "[FullMatrix::SetValue] Row or column number out of range." ) );
        }
        mEigenMatrix ( i,j ) =value;
    }

    T GetValue (int i, int j ) const
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
        {
            throw MathException ( std::string ( "[FullMatrix::GetValue] Row or column number out of range." ) );
        }
        return mEigenMatrix ( i,j );
    }

    void AddValue (int i, int j, const T& value )
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
            throw MathException ( std::string ( "[FullMatrix::AddValue] Row or column number out of range." ) );
        mEigenMatrix ( i,j ) +=value;
    }

    void SubValue (int i, int j, const T& value )
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
            throw MathException ( std::string ( "[FullMatrix::SetValue] Row or column number out of range." ) );
        mEigenMatrix ( i,j )-=value;
    }

    inline int GetNumRows() const
    {
        return mEigenMatrix.rows();
    }

    inline int GetNumColumns() const
    {
        return mEigenMatrix.cols();
    }
#ifndef SWIG
    inline Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& GetEigenMatrix()
    {
        return mEigenMatrix;
    }
    inline const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& GetEigenMatrix() const
    {
        return mEigenMatrix;
    }
#endif
    //! @brief ... convert a sparse matrix into a full matrix of same type
    //! @param sparseMatrix ... sparse matrix
    FullMatrix<T> ( const SparseMatrix<T>& sparseMatrix )
    {
        mEigenMatrix.setZero ( sparseMatrix.GetNumRows(),sparseMatrix.GetNumColumns() );
        sparseMatrix.WriteEntriesToFullMatrix ( *this );
    }

    //! @brief ... copy the matrix
    FullMatrix Copy()
    {
        return ( *this );
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const
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


    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType)
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
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
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
        memcpy ( mEigenMatrix.data(),&(dataVec[0]),numRows * numColumns *sizeof ( T ) );
#ifdef DEBUG_SERIALIZATION
    	std::cout << "finish serialize FullMatrix (load)" << std::endl;
#endif
    }

    //! @brief serializes the class, this is the save routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
#ifdef DEBUG_SERIALIZATION
    	std::cout << "start serialize FullMatrix (save)" << std::endl;
#endif
        std::vector<T> dataVec ( GetNumRows() *GetNumColumns());
        memcpy ( & ( dataVec[0] ),mEigenMatrix.data(),GetNumRows() *GetNumColumns() *sizeof ( T ) );
        int numRows = mEigenMatrix.rows(),
                      numColumns = mEigenMatrix.cols();
        ar & boost::serialization::make_nvp ("Matrix",boost::serialization::base_object< Matrix<T> > ( *this ) )
        & BOOST_SERIALIZATION_NVP(dataVec)
        & BOOST_SERIALIZATION_NVP(numRows)
        & BOOST_SERIALIZATION_NVP(numColumns);
#ifdef DEBUG_SERIALIZATION
    	std::cout << "finish serialize FullMatrix (save)" << std::endl;
#endif
    }


    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
#endif //ENABLE_SERIALIZATION
    //! @brief ... print info about the object
    void Info() const
    {
        Info ( 10,3 );
    }

    //! @brief ... print info about the object
    //! @param width ... total width of the each entry in the matrix when printed
    void Info ( int width ) const
    {
        Info ( width,3 );
    }

    //! @brief ... print info about the object
    //! @param width ... total width of the each entry in the matrix when printed
    //! @param precision ... precision
    void Info ( int width, int precision ) const
    {
        std::cout.precision ( precision );
        for ( int count=0; count<GetNumRows(); count++ )
        {
            for ( int count2=0; count2<GetNumColumns(); count2++ )
            {
                std::cout << std::setw ( width ) << this->Convert2String ( mEigenMatrix.data() [count2*GetNumRows() +count] ) <<" ";
            }
            std::cout<<std::endl;
        }
    }

    //! @brief ... reads a matrix from a file
    //! linesToSkip_ ... number of lines to skip at the beginning
    //! delimiters_ ... delimiters between the entries in each line
    void ReadFromFile ( std::string fileName_, unsigned int linesToSkip_, std::string delimiters_ )
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
                mEigenMatrix.resize ( numRows_,numColumns_ );

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
                        mEigenMatrix ( curRow_,curColumn_ ) =Matrix<T>::ConvertFromString ( *beg );
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
    void ReadFromFile ( std::string fileName_ )
    {
        ReadFromFile ( fileName_,0," " );
    }

    //! @brief ... writes a matrix to a file
    //! delimiter_ ... delimiters between the entries in each line
    void AppendToFile ( std::string fileName_, std::string delimiter_ )
    {
        AppendToFile ( fileName_,delimiter_,"","" );
    }

    //! @brief ... append a matrix to a file
    //! delimiter_ ... delimiters between the entries in each line
    //! stringBefore ... String to write before the matrix (e.g. some comments)
    //! stringAfter ... String to write after the matrix (e.g. some comments)
    void AppendToFile ( std::string fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter )
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
    //! delimiter_ ... delimiters between the entries in each line
    void WriteToFile ( std::string fileName_, std::string delimiter_ )
    {
        WriteToFile ( fileName_,delimiter_,"","" );
    }


    //! @brief ... writes a matrix to a file
    //! delimiter_ ... delimiters between the entries in each line
    //! stringBefore ... String to write before the matrix (e.g. some comments)
    //! stringAfter ... String to write after the matrix (e.g. some comments)
    void WriteToFile ( std::string fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter )
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
    virtual std::string GetTypeId()const;


    //! @brief ... converts any FullMatrix to a DoubleFullMatrix
    //! @return    converted matrix
    FullMatrix<double> Convert2double();

    //! @brief ... converts any FullMatrix to an IntFullMatrix
    //! @return    converted matrix
    FullMatrix<int> Convert2int();


    virtual void ImportFromSLangText ( const char* fileName );

    //! @brief ... imports Matrix from a Vtk ASCII File
    virtual void ImportFromVtkASCIIFile(const char* fileName);

    //! @brief performs a dyadic operator on all matrix entries and the entries of another matrix
    //! @param rOperator        Dyadic Operator
    void Map ( NuTo::DyadicOperator<T>* rDOperator, const FullMatrix<T> &otherMatrix )
    {
        throw MathException ( "FullMatrix::Map - To be implemented." );
    }

    //! @brief performs a monadic operator on all matrix entries
    //! @param rMOperator        Monadic Operator
    virtual void Map ( const NuTo::MonadicOperator<T>* rMOperator )
    {
        for ( int count=0; count<GetNumColumns() *GetNumRows(); count++ )
            mEigenMatrix.data() [count] = rMOperator->Evaluate ( mEigenMatrix.data() [count] );
    }

    //! @brief performs a dyadic operator on all matrix entries with another given value
    //! @param rDOperator        Dyadic Operator
    //! @param rValue ... value
    virtual void Map ( const NuTo::DyadicOperator<T>* rDOperator, const T& rValue )
    {
        for ( int count=0; count<GetNumColumns() *GetNumRows(); count++ )
            mEigenMatrix.data() [count] = rDOperator->Evaluate ( mEigenMatrix.data() [count],rValue );
    }

    //! @brief returns the maximum value of the matrix
    //! if several entries have the same maximum value, only the first one is recovered
    //! @param rRowOutput  row
    //! @param rColumnOutput  column
    virtual T Max ( int& rRowOutput, int& rColumnOutput )
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Max - Maximum for matrix with zero entries cannot be calculated." );

        return mEigenMatrix.maxCoeff ( &rRowOutput,&rColumnOutput );
    }

#ifndef SWIG
//! @brief returns the maximum value of the matrix
    //! if several entries have the same maximum value, only the first one is recovered
    virtual T Max()
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Max - Maximum for matrix with zero entries cannot be calculated." );
        return mEigenMatrix.maxCoeff();
    }
#endif

    //! @brief returns the minimum value of the matrix
    //! if several entries have the same maximum value, only the first one is recovered
    //! @param rRowOutput  row
    //! @param rColumnOutput  column
    virtual T Min ( int& rRowOutput, int& rColumnOutput )
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Min - Minimum for matrix with zero entries cannot be calculated." );

        return mEigenMatrix.minCoeff ( &rRowOutput,&rColumnOutput );
    }

#ifndef SWIG
    //! @brief returns the minimum value of the matrix
    //! if several entries have the same maximum value, only the first one is recovered
    virtual T Min()
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Min - Minimum for matrix with zero entries cannot be calculated." );

        return mEigenMatrix.minCoeff();
    }
#endif

    //! @brief elementwise absolute value of the matrix
    virtual FullMatrix<T> Abs();
//    {
//		  //return FullMatrix<T> ( mEigenMatrix.cwise().abs() );
//    }

    //! @brief elementwise absolute value of the matrix
    virtual FullMatrix<T> Trans()
    {
        return FullMatrix<T> ( mEigenMatrix.transpose() );
    }

    //! @brief Extract a block of the FullMatrix
    FullMatrix<T>  GetBlock ( int rI, int rJ, int rRows, int rCols )
    {
        if ( rI+rRows>mEigenMatrix.rows() )
            throw MathException ( "FullMatrix::GetBlock - rows out of Dimension." );
        if ( rJ+rCols>mEigenMatrix.cols() )
            throw MathException ( "FullMatrix::GetBlock - columns out of Dimension." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::GetBlock - row should not be negative." );
        if ( rJ<0 )
            throw MathException ( "FullMatrix::GetBlock - column should not be negative." );
        return FullMatrix<T> ( mEigenMatrix.block ( rI,rJ,rRows,rCols ) );
    }

    //! @brief Set a block of the FullMatrix
    void  SetBlock ( int rI, int rJ, FullMatrix<T> rBlock )
    {
        if ( rI+rBlock.mEigenMatrix.rows() >mEigenMatrix.rows() )
            throw MathException ( "FullMatrix::SetBlock - rows out of Dimension." );
        if ( rJ+rBlock.mEigenMatrix.cols() >mEigenMatrix.cols() )
            throw MathException ( "FullMatrix::SetBlock - columns out of Dimension." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::SetBlock - row should not be negative." );
        if ( rJ<0 )
            throw MathException ( "FullMatrix::SetBlock - column should not be negative." );
        mEigenMatrix.block ( rI,rJ,rBlock.mEigenMatrix.rows(),rBlock.mEigenMatrix.cols() ) = rBlock.mEigenMatrix;
    }

    //! @brief Extract a row of the FullMatrix
    FullMatrix<T>  GetRow ( int rI )
    {
        if ( rI>=mEigenMatrix.rows() )
            throw MathException ( "FullMatrix::GetRow - row out of Dimension." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::GetRow - row should not be negative." );
        return FullMatrix<T> ( mEigenMatrix.row ( rI ) );
    }

    //! @brief Set a row of the FullMatrix
    void  SetRow ( int rI, FullMatrix<T> rBlock )
    {
        if ( rBlock.mEigenMatrix.rows() !=1 )
            throw MathException ( "FullMatrix::SetRow - Expect a Matrix with a single row as Input." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::SetRow - row should not be negative." );
        if ( rBlock.mEigenMatrix.cols() !=mEigenMatrix.cols() )
            throw MathException ( "FullMatrix::SetRow - number of columns for both matrices must be identical." );
        mEigenMatrix.row ( rI ) = rBlock.mEigenMatrix;
    }

    //! @brief Extract a column of the FullMatrix
    FullMatrix<T>  GetColumn ( int rI )
    {
        if ( rI>=mEigenMatrix.cols() )
            throw MathException ( "FullMatrix::GetColumn - column out of Dimension." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::GetColumn - column should not be negative." );
        return FullMatrix<T> ( mEigenMatrix.col ( rI ) );
    }

    //! @brief Set a column of the FullMatrix
    void  SetColumn ( int rI, FullMatrix<T> rBlock )
    {
        if ( rBlock.mEigenMatrix.cols() !=1 )
            throw MathException ( "FullMatrix::SetColumn - Expect a Matrix with a single column as Input." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::SetColumn - column should not be negative." );
        if ( rBlock.mEigenMatrix.rows() !=mEigenMatrix.rows() )
            throw MathException ( "FullMatrix::SetColumn - number of rows for both matrices must be identical." );
        mEigenMatrix.col ( rI ) = rBlock.mEigenMatrix;
    }

    //! @brief Appends columns to a FullMatrix
    void  AppendColumns ( FullMatrix<T> rBlock )
    {
        if ( rBlock.mEigenMatrix.rows() !=mEigenMatrix.rows() )
            throw MathException ( "FullMatrix::AppendColumns - number of rows for both matrices must be identical." );
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> tmpMatrix ( mEigenMatrix.rows(),mEigenMatrix.cols() +rBlock.mEigenMatrix.cols() );
        tmpMatrix.block ( 0,0,mEigenMatrix.rows(),mEigenMatrix.cols() ) = mEigenMatrix;
        tmpMatrix.block ( 0,mEigenMatrix.cols(),rBlock.mEigenMatrix.rows(),rBlock.mEigenMatrix.cols() ) = rBlock.mEigenMatrix;
        mEigenMatrix = tmpMatrix;
    }

    //! @brief Appends rows to a FullMatrix
    void  AppendRows ( FullMatrix<T> rBlock )
    {
        if ( rBlock.mEigenMatrix.cols() !=mEigenMatrix.cols() )
            throw MathException ( "FullMatrix::AppendRows - number of columns for both matrices must be identical." );
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> tmpMatrix ( mEigenMatrix.rows() +rBlock.mEigenMatrix.rows(),mEigenMatrix.cols() );
        tmpMatrix.block ( 0,0,mEigenMatrix.rows(),mEigenMatrix.cols() ) = mEigenMatrix;
        tmpMatrix.block ( mEigenMatrix.rows(),0,rBlock.mEigenMatrix.rows(),rBlock.mEigenMatrix.cols() ) = rBlock.mEigenMatrix;
        mEigenMatrix = tmpMatrix;
    }

    //! @brief Calculates the norm of a matrix, i.e. for vectors the Euclidean norm
    double  Norm()
    {
        return mEigenMatrix.norm();
    }

    //! @brief Calculates the sum of the entries for each column
    FullMatrix<T> ColumnwiseSum()
    {
        return FullMatrix<T> ( mEigenMatrix.colwise().sum() );
    }

    //! @brief Calculates the sum of the entries for each row
    FullMatrix<T> RowwiseSum()
    {
        return FullMatrix<T> ( mEigenMatrix.rowwise().sum() );
    }

    //! @brief Calculates the minimum of the entries for each column
    FullMatrix<T> ColumnwiseMinCoeff()
    {
        return FullMatrix<T> ( mEigenMatrix.colwise().minCoeff() );
    }

    //! @brief Calculates the minimum of the entries for each row
    FullMatrix<T> RowwiseMinCoeff()
    {
        return FullMatrix<T> ( mEigenMatrix.rowwise().minCoeff() );
    }

    //! @brief Calculates the maximum of the entries for each column
    FullMatrix<T> ColumnwiseMaxCoeff()
    {
        return FullMatrix<T> ( mEigenMatrix.colwise().maxCoeff() );
    }

    //! @brief Calculates the maximum of the entries for each row
    FullMatrix<T> RowwiseMaxCoeff()
    {
        return FullMatrix<T> ( mEigenMatrix.rowwise().maxCoeff() );
    }

    //! @brief Coefficient wise multiplication
    FullMatrix<T> ElementwiseMul ( const FullMatrix<T> &other )
    {
        return FullMatrix<T> ( mEigenMatrix.cwise() *other.mEigenMatrix );
    }

    //! @brief Coefficient wise division
    FullMatrix<T> ElementwiseDiv ( const FullMatrix<T> &other )
    {
        return FullMatrix<T> ( mEigenMatrix.cwise() /other.mEigenMatrix );
    }

    //! @brief Coefficient wise reciprocal
    FullMatrix<T> ElementwiseInverse()
    {
        return FullMatrix<T> ( mEigenMatrix.cwise().inverse() );
    }

#ifndef SWIG
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> mEigenMatrix;
#endif

private:
    //! @brief ... writes a matrix to a file
    //! linesToSkip_ ... number of lines to skip at the beginning
    //! delimiters_ ... delimiters between the entries in each line
    void WriteToFile ( std::ofstream& fileStream, std::string& delimiter_ )
    {
        // go through the matrix and store the values
        for ( int count=0; count<GetNumRows(); count++ )
        {
            for ( int count2=0; count2<GetNumColumns(); count2++ )
            {
                fileStream    <<     this->Convert2String ( mEigenMatrix.data() [count+count2*GetNumRows() ],true,12,18 );
                if ( count2!=GetNumColumns()-1 )
                    fileStream    << delimiter_;
            }
            if ( count<GetNumRows()-1 )
                fileStream    <<  "\n";
        }
    }


};

template<class T>
std::ostream& operator<<(std::ostream& os, const NuTo::FullMatrix<T>& r)
{
    os << r.mEigenMatrix;
return os;
}

} //NAMESPACE NUTO
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::FullMatrix<double>)
BOOST_CLASS_EXPORT_KEY(NuTo::FullMatrix<int>)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif // FULL_MATRIX_H
