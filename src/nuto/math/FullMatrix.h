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

#include "nuto/base/Logger.h"
#include "nuto/math/Matrix.h"
#include "nuto/math/NuToMath.h"
#include "nuto/math/MathException.h"
#include "nuto/math/fortran_routines.h"
#include "nuto/math/SparseMatrix.h"
#include <eigen3/Eigen/Core>


namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date July 2009
//! @brief ... class for full matrices derived from the abstract base class Matrix
template <class T>
class FullMatrix : public Matrix<T>, public Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>

{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif

public:
    //! @brief ... default constructor
    FullMatrix<T>()
    {
        this->resize ( 0,0 );
    }

    //! @brief ... constructor
    //! @param rNumRows ... number of columns
    //! @param rNumColumns ... number of rows
    FullMatrix<T> ( int rNumRows, int rNumColumns)
    {
        if ( rNumRows*rNumColumns>0 )
            this->setZero ( rNumRows,rNumColumns );
        else
            this->resize ( rNumRows,rNumColumns );
    }

#ifndef SWIG
    //! @brief ... constructor
    //! @param rEigenMatrix ... other matrix
    template<typename OtherDerived>
    FullMatrix<T> ( const Eigen::MatrixBase<OtherDerived>& rEigenMatrix): Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>(rEigenMatrix)
    {
    }
#endif

    //! @brief ... constructur
    //! @param numRows_ ... number of columns
    //! @param numColumns_ ... number of rows
    //! @param entries_ ... vector containing the matrix in column-major orientation
    FullMatrix<T> ( int numRows_, int numColumns_, const std::vector<T>& entries_ )
    {
        this->resize ( numRows_,numColumns_ );
        const T *ptr = &entries_[0];
        for ( int j=0; j<this->cols(); ++j )              // loop over columns
            for ( int i=0; i<this->rows(); ++i, ptr++ )   // loop over rows
                (*this)( i,j ) = *ptr;                           // to access matrix coefficients,
    }

    //! @brief ... constructor
    //! @brief ... creates a FullVector(means a FullMatrix size(m,1) )
    //! @param entries_ ... vector containing the matrix in column-major orientation
    FullMatrix<T> ( const std::vector<T>& entries_ )
    {
        this->resize ( entries_.size(),1 );
        const T *ptr = &entries_[0];

		for ( int i=0; i<this->rows(); i++, ptr++ )   	// loop over rows
			(*this)( i,0 ) = *ptr;                     	// to access matrix coefficients,
    }

   //! @brief ... copy constructor
    //! @param  rOther ... copied element
    FullMatrix<T> ( const FullMatrix<T>& rOther ): Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>(rOther)
    {
    }

    //! @brief ... assignment constructor
	//! @param  rOther ... copied element
	FullMatrix<T>& operator=( const FullMatrix<T>& rOther )
	{
		if (this != &rOther)
		{
			this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator=(rOther);
		}
		return *this;
	}

	//! @brief ... assignment constructor
	//! @param  rOther ... copied element
    template<typename OtherDerived>
	FullMatrix<T>& operator=( const Eigen::MatrixBase <OtherDerived>& other)
	{
    	this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator=(other);
    	return *this;
	}

    //! @brief ... resize matrix (everything is deleted)
    //! @param rows ... number of rows
    //! @param cols ... number of columns
    void Resize ( int rows, int cols )
    {
    	if ( rows*cols>0 )
    		this->setZero ( rows, cols );
        else
        	this->resize ( rows,cols );
    }

    //! @brief ... resize matrix (nothing is deleted)
    //! @param rows ... number of rows
    //! @param cols ... number of columns
    void ConservativeResize ( int rows, int cols )
    {
        this->conservativeResize ( rows,cols );
    }

    //! @brief ... resize matrix (nothing is deleted)
    //! @param rows ... number of rows
    void ConservativeResizeRows (int rows )
    {
    	this->conservativeResize(rows, Eigen::NoChange_t());
    }

    //! @brief ... resize matrix (nothing is deleted)
    //! @param rows ... number of columns
    void ConservativeResizeCols (int cols )
    {
    	this->conservativeResize(Eigen::NoChange_t(), cols);
    }


    //! @brief ... add another matrix to this matrix and return the result matrix
    //! @param other ... other matrix
    //! @return resultMatrix = thisMatrix + otherMatrix
    FullMatrix<T> operator+ ( const FullMatrix<T> &other ) const
    {
        if ( this->rows() !=other.rows() || this->cols() !=other.cols() )
            throw MathException ( std::string ( "[FullMatrix::operator+] Row or column number must be identical." ) );
        return this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator+(other);
    }

    //! @brief ... subtract another matrix from this matrix and return result matrix
    //! @param other ... other matrix
    //! @return resultMatrix = thisMatrix - otherMatrix
    FullMatrix<T> operator- ( const FullMatrix<T> &other ) const
    {
        if ( this->rows() !=other.rows() || this->cols() !=other.cols() )
            throw MathException ( std::string ( "[FullMatrix::operator+] Row or column number must be identical." ) );
        return this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator-(other);
    }

    //! @brief ... multiply this matrix by another matrix and return the result matrix
    //! @param other ... other matrix
    //! @return resultMatrix = thisMatrix * otherMatrix
    FullMatrix<T> operator* ( const FullMatrix<T> &other ) const
    {
        if ( this->cols() !=other.rows() )
            throw MathException ( std::string ( "[FullMatrix::operator*] Number of columns of the first matrix must be identical to the number of rows of the second matrix." ) );
        return this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator*(other);
    }

    //! @brief ... multiply this matrix by a vector (which is stored in a vector
    //! @param other ... vector
    //! @return matrix
    FullMatrix<T> operator* ( const std::vector<T> &other ) const;

    //! @brief ... multiply this matrix by a scalar factor and return the result matrix
    //! @param other ... scalar factor
    //! @return resultMatrix = scalarFactor * thisMatrix
    FullMatrix<T> operator* ( const T &other ) const
    {
    	return this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator*(other);
    }

    //! @brief ... add in place another matrix to this matrix and return a reference to this matrix
    //! @param other ... other matrix
    //! @return reference to this matrix
    FullMatrix<T>& operator+= ( const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &other )
    {
        if ( this->rows() !=other.rows() || this->cols() !=other.cols() )
            throw MathException ( std::string ( "[FullMatrix::operator+=] Row or column number must be identical." ) );
        this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator+=(other);
        return *this;
    }

    //! @brief ... subtract in place another matrix from this matrix and return reference to this matrix
    //! @param other ... other matrix
    //! @return reference to this matrix
    FullMatrix<T>& operator-= ( const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> &other )
	{
		if ( this->rows() !=other.rows() || this->cols() !=other.cols() )
			throw MathException ( std::string ( "[FullMatrix::operator+=] Row or column number must be identical." ) );
		this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator-=(other);
        return *this;
	}
    //! @brief ... scale this matrix in place by a scalar factor and return a reference to this matrix
    //! @param other ... scalar factor
    //! @return reference to this matrix
    FullMatrix<T>& operator*= ( const T &other)
    {
    	this->Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::operator*=(other);
        return *this;
    }

    //! @brief ... add a value to all entries of the matrix and return a reference to this matrix
    //! @param other ... value to be added
    //! @return reference to this matrix
    FullMatrix<T>& operator+= ( const T &other)
    {
        this->array()+=other;
        return *this;
    }
/*

    //! @brief ... access operator
    //! @param i ... row
    //! @param j ... column
    //! @return reference to the (i,j)-th entry of the matrix
    T& operator() ( int i, int j )
    {
        if ( i>=mEigenMatrix.rows() || j>=mEigenMatrix.cols() || i<0 || j<0 )
            throw MathException ( std::string ( "[FullMatrix::operator()] Row or column number out of range." ) );
        return  mEigenMatrix ( i,j );
    }

#ifndef SWIG
    //! @brief ... access operator
    //! @param i ... row
    //! @param j ... column
    //! @return the (i,j)-th entry of the matrix
    const T operator() ( int i, int j ) const
    {
        if ( i>=mEigenMatrix.rows() || j>=mEigenMatrix.cols() || i<0 || j<0 )
            throw MathException ( std::string ( "[FullMatrix::operator()] Row or column number out of range." ) );
        return  mEigenMatrix ( i,j );
    }
#endif
*/
    //! @brief ... calculates the scalar product
    //! @param other ... other matrix
    //! @return reference to this matrix
    T Dot( const FullMatrix<T> &other )
    {
        if ( this->rows() !=other.rows() || this->cols() !=1 || other.cols()!=1 )
            throw MathException ( std::string ( "[FullMatrix::Dot] Dot product requires number of rows to be identical (and only a single column)." ) );
         return (this->transpose()*other).eval()(0,0);
    }

    //! @brief ... set the value of the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @param value ... value of the matrix entry
    void SetValue (int i, int j, const T& value )
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
        {
            throw MathException ( std::string ( "[FullMatrix::SetValue] Row or column number out of range." ) );
        }
        (*this)( i,j ) =value;
    }

    //! @brief ... get the value of the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @return the value of the (i,j)-th matrix entry
    T GetValue (int i, int j ) const
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
        {
            throw MathException ( std::string ( "[FullMatrix::GetValue] Row or column number out of range." ) );
        }
        return (*this) ( i,j );
    }

    //! @brief ... add a value to the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @param value ... value which is added to the (i,j)-th matrix entry
    void AddValue (int i, int j, const T& value )
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
            throw MathException ( std::string ( "[FullMatrix::AddValue] Row or column number out of range." ) );
        (*this) ( i,j ) +=value;
    }

    //! @brief ... subtract a value from the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @param value ... value which is subtracted from the (i,j)-th matrix entry
    void SubValue (int i, int j, const T& value )
    {
        if ( i>=GetNumRows() || i<0 || j>=GetNumColumns() || j<0 )
            throw MathException ( std::string ( "[FullMatrix::SetValue] Row or column number out of range." ) );
        (*this) ( i,j )-=value;
    }

    //! @brief ... get number of rows
    //! @return number of rows
    inline int GetNumRows() const
    {
        return this->rows();
    }

    //! @brief ... get number of columns
    //! @return number of columns
    inline int GetNumColumns() const
    {
        return this->cols();
    }

/*
#ifndef SWIG
    //! @brief get a reference to the eigen matrix
    //! @return reference to eigen matrix
    inline Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& GetEigenMatrix()
    {
        return mEigenMatrix;
    }

    //! @brief get a reference to the eigen matrix
    //! @return reference to eigen matrix
    inline const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& GetEigenMatrix() const
    {
        return mEigenMatrix;
    }
#endif
*/
    //! @brief ... convert a sparse matrix into a full matrix of same type
    //! @param sparseMatrix ... sparse matrix
    FullMatrix<T> ( const SparseMatrix<T>& sparseMatrix )
    {
    	this->setZero ( sparseMatrix.GetNumRows(),sparseMatrix.GetNumColumns() );
        sparseMatrix.WriteEntriesToFullMatrix ( *this );
    }

    //! @brief ... copy the matrix
    //! @return a copy of this matrix
    FullMatrix Copy()
    {
        return ( *this );
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType) const;

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType);

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
        memcpy ( this->data(),&(dataVec[0]),numRows * numColumns *sizeof ( T ) );
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


    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
#endif //ENABLE_SERIALIZATION
    //! @brief ... print info about the object
    void Info() const
    {
        Info ( 10,3 );
    }

    //! @brief ... print info about the object
    //! @param width ... total width of each entry in the matrix when printed
    void Info ( int width ) const
    {
        Info ( width,3 );
    }

    //! @brief ... print info about the object
    //! @param width ... total width of the each entry in the matrix when printed
    //! @param precision ... precision (number of fractional digits)
    //! @param rScientific ... switch for scientific notation
    void Info ( int rWidth, int rPrecision, bool rScientific = false ) const
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
    void Out(NuTo::Logger& rLogger, int rWidth, int rPrecision, bool rScientific = false )const
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
    void ReadFromFile ( std::string fileName_ )
    {
        ReadFromFile ( fileName_,0," " );
    }

    //! @brief ... writes a matrix to a file
    //! @param fileName_ ... file name
    //! @param delimiter_ ... delimiters between the entries in each line
    void AppendToFile ( std::string fileName_, std::string delimiter_ ) const
    {
        AppendToFile ( fileName_,delimiter_,"","" );
    }

    //! @brief ... append a matrix to a file
    //! @param fileName_ ... file name
    //! @param delimiter_ ... delimiters between the entries in each line
    //! @param stringBefore ... String to write before the matrix (e.g. some comments)
    //! @param stringAfter ... String to write after the matrix (e.g. some comments)
    void AppendToFile ( std::string fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter ) const
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
    void WriteToFile ( std::string fileName_, std::string delimiter_ ) const
    {
        WriteToFile ( fileName_,delimiter_,"","" );
    }


    //! @brief ... writes a matrix to a file
    //! @param fileName_ ... file name
    //! @param delimiter_ ... delimiters between the entries in each line
    //! @param stringBefore ... String to write before the matrix (e.g. some comments)
    //! @param stringAfter ... String to write after the matrix (e.g. some comments)
    void WriteToFile ( std::string fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter ) const
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

    //! @brief solves a system of linear equations with a symmetric positive definite coefficient matrix using Cholesky factorization (LAPACK)
    //! @param rRHS ... matrix of right-hand-side vectors
    //! @param rSolution ... solution matrix
    void SolveCholeskyLapack(const FullMatrix<double>& rRHS, FullMatrix<double>& rSolution) const;

    //! @brief calculates the inverse of a symmetric positive definite matrix using Cholesky factorization (LAPACK)
    //! @param rInverse ... inverse matrix
    void InverseCholeskyLapack(FullMatrix<double>& rInverse) const;

    //! @brief calculates the inverse of a matrix
    //! @param rInverse ... inverse matrix
    FullMatrix<T> Inverse() const;

    //! @brief calculates the eigenvalues
    //! @param rEigenValues ... eigenvalues
    void EigenValuesSymmetric(FullMatrix<double>& rEigenValues)const;

    //! @brief calculates the eigenvectors
    //! @param rEigenVectors ... eigenvectors
    void EigenVectorsSymmetric(FullMatrix<double>& rEigenVectors)const;

    //! @brief ... imports a matrix from a SLang ASCII file
    //! @param fileName ... file name
    virtual void ImportFromSLangText ( const char* fileName );

    //! @brief ... imports Matrix from a Vtk ASCII File
    //! @param fileName ... file name
    virtual void ImportFromVtkASCIIFile(const char* fileName);

    //! @brief performs a dyadic operator on all matrix entries and the entries of another matrix
    //! @param rDOperator ... Dyadic Operator
    //! @param otherMatrix ... other matrix
    void Map ( NuTo::DyadicOperator<T>* rDOperator, const FullMatrix<T> &otherMatrix )
    {
        throw MathException ( "FullMatrix::Map - To be implemented." );
    }

    //! @brief performs a monadic operator on all matrix entries
    //! @param rMOperator ... Monadic Operator
    virtual void Map ( const NuTo::MonadicOperator<T>* rMOperator )
    {
        for ( int count=0; count<GetNumColumns() *GetNumRows(); count++ )
        	this->data() [count] = rMOperator->Evaluate ( this->data() [count] );
    }

    //! @brief performs a dyadic operator on all matrix entries with another given value
    //! @param rDOperator ... Dyadic Operator
    //! @param rValue ... value
    virtual void Map ( const NuTo::DyadicOperator<T>* rDOperator, const T& rValue )
    {
        for ( int count=0; count<GetNumColumns() *GetNumRows(); count++ )
            this->data() [count] = rDOperator->Evaluate ( this->data() [count],rValue );
    }

    //! @brief returns the maximum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @param rRowOutput ... row
    //! @param rColumnOutput ... column
    //! @return maximum value of the matrix
    virtual T Max ( int& rRowOutput, int& rColumnOutput ) const
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Max - Maximum for matrix with zero entries cannot be calculated." );

        return this->maxCoeff ( &rRowOutput,&rColumnOutput );
    }

#ifndef SWIG
    //! @brief returns the maximum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @return maximum value of the matrix
    virtual T Max() const
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Max - Maximum for matrix with zero entries cannot be calculated." );
        return this->maxCoeff();
    }
#endif

    //! @brief returns the minimum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @param rRowOutput ... row
    //! @param rColumnOutput ...  column
    //! @return minimum value of the matrix
    virtual T Min ( int& rRowOutput, int& rColumnOutput ) const
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Min - Minimum for matrix with zero entries cannot be calculated." );

        return this->minCoeff ( &rRowOutput,&rColumnOutput );
    }

#ifndef SWIG
    //! @brief returns the minimum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @return minimum value of the matrix
    virtual T Min() const
    {
        if ( GetNumColumns() ==0 || GetNumRows() ==0 )
            throw MathException ( "FullMatrix::Min - Minimum for matrix with zero entries cannot be calculated." );

        return this->minCoeff();
    }
#endif

    //! @brief elementwise absolute value of the matrix
    virtual FullMatrix<T> Abs() const;

    //! @brief calcualte the transpose of this matrix
    //! @return transpose of the matrix
    virtual FullMatrix<T> Trans() const
    {
        return FullMatrix<T> ( this->transpose() );
    }

    //! @brief ... extract a block (submatrix) from this matrix
    //! @param rI ... start row
    //! @param rJ ... start column
    //! @param rRows ... number of rows of the submatrix
    //! @param rCols ... number of columns of the submatrix
    //! @return a submatrix with size (rRows,rCols) extracted at (rI,rJ)-entry of this matrix
    FullMatrix<T>  GetBlock ( int rI, int rJ, int rRows, int rCols ) const
    {
        if ( rI+rRows>this->rows() )
            throw MathException ( "FullMatrix::GetBlock - rows out of Dimension." );
        if ( rJ+rCols>this->cols() )
            throw MathException ( "FullMatrix::GetBlock - columns out of Dimension." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::GetBlock - row should not be negative." );
        if ( rJ<0 )
            throw MathException ( "FullMatrix::GetBlock - column should not be negative." );
        return FullMatrix<T> ( this->block ( rI,rJ,rRows,rCols ) );
    }

    //! @brief ... set a block (submatrix) in this Matrix
    //! @param rI ... start row
    //! @param rJ ... start column
    //! @param rBlock ... submatrix
    void  SetBlock ( int rI, int rJ, const FullMatrix<T>& rBlock )
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
    void  AddBlock ( int rI, int rJ, const FullMatrix<T>& rBlock )
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
    FullMatrix<T>  GetRow ( int rI ) const
    {
        if ( rI>=this->rows() )
            throw MathException ( "FullMatrix::GetRow - row out of Dimension." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::GetRow - row should not be negative." );
        return FullMatrix<T> ( this->row ( rI ) );
    }

    //! @brief ... set a row in this FullMatrix
    //! @param rI ... row
    //! @param rBlock ... new row intries
    void  SetRow ( int rI, const FullMatrix<T>& rBlock )
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
    FullMatrix<T>  GetColumn ( int rI ) const
    {
        if ( rI>=this->cols() )
            throw MathException ( "FullMatrix::GetColumn - column out of Dimension." );
        if ( rI<0 )
            throw MathException ( "FullMatrix::GetColumn - column should not be negative." );
        return FullMatrix<T> ( this->col ( rI ) );
    }

    //! @brief ... set a column in this matrix
    //! @param rI ... column
    //! @param rBlock ... new column entries
    void  SetColumn ( int rI, const FullMatrix<T>& rBlock )
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
    void  AppendColumns ( FullMatrix<T> rBlock )
    {
        if ( rBlock.rows() !=this->rows() )
            throw MathException ( "FullMatrix::AppendColumns - number of rows for both matrices must be identical." );
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> tmpMatrix ( this->rows(),this->cols() +rBlock.cols() );
        tmpMatrix.block ( 0,0,this->rows(),this->cols() ) = *this;
        tmpMatrix.block ( 0,this->cols(),rBlock.rows(),rBlock.cols() ) = rBlock;
        *this = tmpMatrix;
    }

    //! @brief ... appends rows to this matrix
    //! @param rBlock ... matrix storing the rows
    void  AppendRows ( FullMatrix<T> rBlock )
    {
        if ( rBlock.cols() !=this->cols() )
            throw MathException ( "FullMatrix::AppendRows - number of columns for both matrices must be identical." );
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> tmpMatrix ( this->rows() +rBlock.rows(),this->cols() );
        tmpMatrix.block ( 0,0,this->rows(),this->cols() ) = (*this);
        tmpMatrix.block ( this->rows(),0,rBlock.rows(),rBlock.cols() ) = rBlock;
        *this = tmpMatrix;
    }

    //! @brief ... calculates the norm of this matrix, i.e. for vectors the Euclidean norm
    //! @return norm of this matrix
    double  Norm() const;

    //! @brief ... calculates the sum of the entries for each column
    //! @return sum of the entries for each column
    FullMatrix<T> ColumnwiseSum() const
    {
        return FullMatrix<T> ( this->colwise().sum() );
    }

    //! @brief ... calculates the sum of the entries for each row
    //! @return sum of the entries for each row
    FullMatrix<T> RowwiseSum() const
    {
        return FullMatrix<T> ( this->rowwise().sum() );
    }

    //! @brief ... calculates the minimum of the entries for each column
    //! @return the minimum of the entries for each column
    FullMatrix<T> ColumnwiseMinCoeff() const
    {
        return FullMatrix<T> ( this->colwise().minCoeff() );
    }

    //! @brief ... calculates the minimum of the entries for each row
    //! @return the minimum of the entries for each row
    FullMatrix<T> RowwiseMinCoeff() const
    {
        return FullMatrix<T> ( this->rowwise().minCoeff() );
    }

    //! @brief ... calculates the maximum of the entries for each column
    //! @return the maximum of the entries for each column
    FullMatrix<T> ColumnwiseMaxCoeff() const
    {
        return FullMatrix<T> ( this->colwise().maxCoeff() );
    }

    //! @brief ... calculates the maximum of the entries for each row
    //! @return the maximum of the entries for each row
    FullMatrix<T> RowwiseMaxCoeff() const
    {
        return FullMatrix<T> ( this->rowwise().maxCoeff() );
    }

    //! @brief ... coefficient wise multiplication of this matrix with another matrix
    //! @param other ... other matrix
    //! @return a matrix which is obtained by a coefficient wise multiplication of this matrix with another matrix
    FullMatrix<T> ElementwiseMul ( const FullMatrix<T> &other ) const
    {
        return FullMatrix<T> ( (this->array() *other.array()).matrix() );
    }

    //! @brief ... coefficient wise division of this matrix by another matrix
    //! @param other ... other matrix
    //! @return a matrix which is obtained by a coefficient wise division of this matrix by another matrix
    FullMatrix<T> ElementwiseDiv ( const FullMatrix<T> &other ) const
    {
        return FullMatrix<T> ( (this->array() /other.array()).matrix() );
    }

    //! @brief ... coefficient wise reciprocal
    //! @return a matrix which is obtained by coefficient wise reciprocal of this matrix
    FullMatrix<T> ElementwiseInverse() const
    {
        return FullMatrix<T> ( this->array().inverse().matrix());
    }

    //! @brief ... sorts the rows of a matrix based on the entries in column rCol
    FullMatrix<T> SortRow(int rCol)
    {
    	if (rCol<0 || rCol>=this->GetNumColumns())
    		throw MathException("[FullMatrix<T>::SortRow] incorrect row");

    	struct mySortClassRow {
    	    const FullMatrix<T> *mat;
    	    int colnum;
    	    mySortClassRow(const FullMatrix<T> &m, int keycol)  : mat(&m), colnum(keycol){};
    	    bool operator()(const int &lhs, const int &rhs) const
    	    {
    	    	return (*mat)(lhs,colnum) < (*mat)(rhs,colnum);
    	    }
    	};

    	std::vector<unsigned int> rows(this->GetNumRows());
    	for (unsigned int count=0; count<rows.size(); count++)
    	{
    		rows[count]=count;
    	}

    	//sort the index vector according to the entries in the matrix
    	std::sort(rows.begin(),rows.end(),mySortClassRow(*this,rCol));

        //rearrange the full matrix rows
    	//for simplicity, just copy the matrix
    	NuTo::FullMatrix<T> myOrigData(*this);
    	for (unsigned int count=0; count<rows.size(); count++)
    	{
    		for ( int count2=0; count2<this->GetNumColumns(); count2++)
    		{
    			(*this)(count,count2) = myOrigData(rows[count],count2);
    		}
    	}
    	return *this;
    }

//#ifndef SWIG
//    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> mEigenMatrix; //!< matrix data are stored as eigen matrix
//#endif

private:
    //! @brief ... writes a matrix to a file
    //! @param fileStream ... output file stream
    //! @param delimiter_ ... delimiters between the entries in each line
    void WriteToFile ( std::ofstream& fileStream, std::string& delimiter_ ) const
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


};

template<class T>
std::ostream& operator<<(std::ostream& os, const NuTo::FullMatrix<T>& r)
{
    os << (Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>&) r;
    return os;
}

} //NAMESPACE NUTO
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::FullMatrix<double>)
BOOST_CLASS_EXPORT_KEY(NuTo::FullMatrix<int>)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

//serilization function for eigen matrices
template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void serialize(
	Archive & ar,
	Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
	const unsigned int file_version
)
{
	ar & boost::serialization::make_array(t.data(), t.size());
}

#endif // FULL_MATRIX_H
