// $Id: FullMatrix.h 623 2013-04-05 08:14:22Z unger3 $

#ifndef FULL_MATRIX_DEF_H
#define FULL_MATRIX_DEF_H

//plugins for the eigen matrices
#define EIGEN_MATRIXBASE_PLUGIN "MatrixBaseAddons.h"

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/array.hpp>
#include <boost/utility/identity_type.hpp>
#endif  // ENABLE_SERIALIZATION


#include "nuto/math/Matrix.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{
//forward declarations
template<typename T> std::string GetTypeIdBaseDataType();
template <class T> class SparseMatrix;
template <class T,int rows> class FullVector;
class Logger;

//! @author Jörg F. Unger, ISM
//! @date July 2009
//! @brief ... class for full matrices derived from the abstract base class Matrix
template <class T, int rows=Eigen::Dynamic, int cols=Eigen::Dynamic>
class FullMatrix : public Matrix<T>, public Eigen::Matrix<T, rows, cols>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif

public:
    //! @brief ... default constructor
    FullMatrix<T,rows,cols>();

    //! @brief ... constructor
    //! @param rNumRows ... number of columns
    //! @param rNumColumns ... number of rows
    FullMatrix<T,rows,cols> ( int rNumRows, int rNumColumns);

#ifdef ENABLE_NUMPY
    void convrtNumpyToMatrix(T * inData,int rRow,int rCol);
    void convrtMatrixToNumpy(T * indata,int rRow,int rCol);
#endif
     
#ifndef SWIG
    //! @brief ... constructor
    //! @param rEigenMatrix ... other matrix
    template<typename OtherDerived>
    FullMatrix<T,rows,cols> ( const Eigen::MatrixBase<OtherDerived>& rOther);
#endif

    //! @brief ... constructur
    //! @param numRows_ ... number of columns
    //! @param numColumns_ ... number of rows
    //! @param entries_ ... vector containing the matrix in column-major orientation
    FullMatrix<T,rows,cols> ( int numRows_, int numColumns_, const std::vector<T>& entries_ );

   //! @brief ... copy constructor
    //! @param  rOther ... copied element
    FullMatrix<T,rows,cols> ( const FullMatrix<T,rows,cols>& rOther );

    //! @brief ... copy constructor
    //! @param  rOther ... copied element
    FullMatrix<T,rows,cols> ( const SparseMatrix<T>& rOther );

#ifndef SWIG
    //! @brief ... assignment constructor
	//! @param  rOther ... copied element
    template<typename OtherDerived>
	FullMatrix<T,rows,cols>& operator=( const Eigen::MatrixBase <OtherDerived>& other);


    //! @brief ... add in place another matrix to this matrix and return a reference to this matrix
    //! @param other ... other matrix
    //! @return reference to this matrix
    template<typename OtherDerived>
	FullMatrix<T,rows,cols>& operator+=( const Eigen::MatrixBase <OtherDerived>& other);

    //! @brief ... subtract in place another matrix from this matrix and return reference to this matrix
    //! @param other ... other matrix
    //! @return reference to this matrix
    template<typename OtherDerived>
	FullMatrix<T,rows,cols>& operator-=( const Eigen::MatrixBase <OtherDerived>& other);

    //! @brief ... assignment constructor for diagonal matrices
    //! @param  rOther ... copied element
    template<typename OtherDerived>
    FullMatrix<T,rows,cols>& operator=( const Eigen::DiagonalBase< OtherDerived > &other);

    //! @brief ... assignment constructor for diagonal matrices
    //! @param  rOther ... copied element
    template<typename OtherDerived>
    FullMatrix<T,rows,cols>& operator+=( const Eigen::DiagonalBase< OtherDerived > &other);

#endif

    //! @brief ... resize matrix (everything is deleted)
    //! @param rows ... number of rows
    //! @param cols ... number of columns
    void Resize ( int rRows, int rCols ) override;

#ifndef SWIG
    //! @brief ... resize matrix (nothing is deleted)
    //! @param rows ... number of rows
    //! @param cols ... number of columns
    void ConservativeResize ( int rRows, int rCols );

    //! @brief ... resize matrix (nothing is deleted)
    //! @param rows ... number of rows
    void ConservativeResizeRows (int rRows );

    //! @brief ... resize matrix (nothing is deleted)
    //! @param rows ... number of columns
    void ConservativeResizeCols (int rCols );
#endif

    //! @brief ... scale this matrix in place by a scalar factor and return a reference to this matrix
    //! @param other ... scalar factor
    //! @return reference to this matrix
    FullMatrix<T,rows,cols>& operator*= ( const T &other);

    //! @brief ... add a value to all entries of the matrix and return a reference to this matrix
    //! @param other ... value to be added
    //! @return reference to this matrix
    FullMatrix<T,rows,cols>& operator+= ( const T &other);

    //! @brief ... calculates the scalar product
    //! @param other ... other matrix
    //! @return reference to this matrix
    T Dot( const FullMatrix<T,rows,cols> &other );

    //! @brief ... set the value of the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @param value ... value of the matrix entry
    void SetValue (int i, int j, const T& value );

    //! @brief ... set the value of the (i)-th vector entry
    //! @param i ... row
    //! @param value ... value of the matrix entry
    void SetValue (int i, const T& value );

    //! @brief ... get the value of the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @return the value of the (i,j)-th matrix entry
    T GetValue (int i, int j) const;

    //! @brief ... get the value of the (i)-th vector entry
    //! @param i ... row
    //! @return the value of the (i,j)-th matrix entry
    T GetValue (int i) const;

    //! @brief ... add a value to the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @param value ... value which is added to the (i,j)-th matrix entry
    void AddValue (int i, int j, const T& value );

    //! @brief ... subtract a value from the (i,j)-th matrix entry
    //! @param i ... row
    //! @param j ... column
    //! @param value ... value which is subtracted from the (i,j)-th matrix entry
    void SubValue (int i, int j, const T& value );

    //! @brief ... get number of rows
    //! @return number of rows
    int GetNumRows() const;

    //! @brief ... get number of columns
    //! @return number of columns
    int GetNumColumns() const;

#ifdef ENABLE_SERIALIZATION
    void Save ( const std::string &filename, std::string rType)const;

    void Restore ( const std::string &filename,  std::string rType);
#ifndef SWIG

    //! @brief serializes the class, this is the load routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    //! @brief serializes the class, this is the save routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

     BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
#endif //ENABLE_SERIALIZATION

    //! @brief ... print info about the object
    void Info() const;

    //! @brief ... print info about the object
    //! @param width ... total width of each entry in the matrix when printed
    void Info ( int width ) const;

    //! @brief ... print info about the object
    //! @param width ... total width of the each entry in the matrix when printed
    //! @param precision ... precision (number of fractional digits)
    //! @param rScientific ... switch for scientific notation
    void Info ( int rWidth, int rPrecision, bool rScientific = false ) const;

    //! @brief ... print info about the object into logger class
    //! @param width ... total width of the each entry in the matrix when printed
    //! @param precision ... precision (number of fractional digits)
    //! @param rScientific ... switch for scientific notation
    void Out(NuTo::Logger& rLogger, int rWidth, int rPrecision, bool rScientific = false )const;

    //! @brief ... reads a matrix from a file
    //! fileName_ ... file name
    //! linesToSkip_ ... number of lines to skip at the beginning
    //! delimiters_ ... delimiters between the entries in each line
    void ReadFromFile ( std::string fileName_, unsigned int linesToSkip_, std::string delimiters_ );

    //! @brief ... reads a matrix from a file
    //! @param  fileName_ ... file name
    void ReadFromFile ( std::string fileName_ );

    //! @brief ... writes a matrix to a file
    //! @param fileName_ ... file name
    //! @param delimiter_ ... delimiters between the entries in each line
    void AppendToFile ( const std::string& fileName_, std::string delimiter_ ) const;

    //! @brief ... append a matrix to a file
    //! @param fileName_ ... file name
    //! @param delimiter_ ... delimiters between the entries in each line
    //! @param stringBefore ... String to write before the matrix (e.g. some comments)
    //! @param stringAfter ... String to write after the matrix (e.g. some comments)
    void AppendToFile ( const std::string& fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter ) const;

    //! @brief ... writes a matrix to a file
    //! @param fileName_ ... file name
    //! @param delimiter_ ... delimiters between the entries in each line
    void WriteToFile ( const std::string& fileName_, std::string delimiter_ ) const;


    //! @brief ... writes a matrix to a file
    //! @param fileName_ ... file name
    //! @param delimiter_ ... delimiters between the entries in each line
    //! @param stringBefore ... String to write before the matrix (e.g. some comments)
    //! @param stringAfter ... String to write after the matrix (e.g. some comments)
    void WriteToFile ( const std::string& fileName_, std::string delimiter_, std::string stringBefore, std::string stringAfter ) const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const;

    //! @brief ... converts any FullMatrix to a DoubleFullMatrix
    //! @return    converted matrix
    FullMatrix<double,rows,cols> Convert2double();

    //! @brief ... converts any FullMatrix to a IntFullMatrix
    //! @return    converted matrix
    FullMatrix<int,rows,cols> Convert2int();

    //! @brief solves a system of linear equations with a symmetric positive definite coefficient matrix using Cholesky factorization (LAPACK)
    //! @param rRHS ... matrix of right-hand-side vectors
    //! @param rSolution ... solution matrix
    void SolveCholeskyLapack(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRHS, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSolution) const;

    //! @brief calculates the inverse of a symmetric positive definite matrix using Cholesky factorization (LAPACK)
    //! @param rInverse ... inverse matrix
    void InverseCholeskyLapack(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInverse) const;

    //! @brief calculates the inverse of a matrix
    //! @param rInverse ... inverse matrix
    FullMatrix<T,rows,cols> Inverse() const;

    //! @brief calculates the eigenvalues
    //! @param rEigenValues ... eigenvalues
    void EigenValuesSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues) const;

    //! @brief calculates the eigenvectors
    //! @param rEigenVectors ... eigenvectors
    void EigenVectorsSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors)const;

    //! @brief calculates the eigenvectors
    //! @param rEigenVectors ... eigenvectors
    void GeneralizedEigenVectorsSymmetric(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rM, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors)const;

    //! @brief ... imports a matrix from a SLang ASCII file
    //! @param fileName ... file name
    //void ImportFromSLangText ( const char* fileName );

    //! @brief ... imports Matrix from a Vtk ASCII File
    //! @param fileName ... file name
    void ImportFromVtkASCIIFile(const char* fileName);

    //! @brief performs a dyadic operator on all matrix entries and the entries of another matrix
    //! @param rDOperator ... Dyadic Operator
    //! @param otherMatrix ... other matrix
    void Map ( NuTo::DyadicOperator<T>* rDOperator, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &otherMatrix );

    //! @brief performs a monadic operator on all matrix entries
    //! @param rMOperator ... Monadic Operator
    void Map ( const NuTo::MonadicOperator<T>* rMOperator );

    //! @brief performs a dyadic operator on all matrix entries with another given value
    //! @param rDOperator ... Dyadic Operator
    //! @param rValue ... value
    void Map ( const NuTo::DyadicOperator<T>* rDOperator, const T& rValue );

#ifndef SWIG
    //! @brief returns the maximum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @param rRowOutput ... row
    //! @param rColumnOutput ... column
    //! @return maximum value of the matrix
    T Max ( int& rRowOutput, int& rColumnOutput ) const;
#endif

    //! @brief returns the maximum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @return maximum value of the matrix
    T Max() const;

#ifndef SWIG
    //! @brief returns the minimum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @param rRowOutput ... row
    //! @param rColumnOutput ...  column
    //! @return minimum value of the matrix
    T Min ( int& rRowOutput, int& rColumnOutput ) const;
#endif

    //! @brief returns the minimum value of the matrix (if several entries have the same maximum value, only the first one is recovered)
    //! @return minimum value of the matrix
    T Min() const;

    //! @brief elementwise absolute value of the matrix
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> Abs() const;

    //! @brief calcualte the transpose of this matrix
    //! @return transpose of the matrix
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> Trans() const;

    //! @brief ... extract a block (submatrix) from this matrix
    //! @param rI ... start row
    //! @param rJ ... start column
    //! @param rRows ... number of rows of the submatrix
    //! @param rCols ... number of columns of the submatrix
    //! @return a submatrix with size (rRows,rCols) extracted at (rI,rJ)-entry of this matrix
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>  GetBlock ( int rI, int rJ, int rRows, int rCols ) const;

    //! @brief ... set a block (submatrix) in this Matrix
    //! @param rI ... start row
    //! @param rJ ... start column
    //! @param rBlock ... submatrix
    void  SetBlock ( int rI, int rJ, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock );

    //! @brief ... add a block (submatrix) in this Matrix
    //! @param rI ... start row
    //! @param rJ ... start column
    //! @param rBlock ... submatrix
    void  AddBlock ( int rI, int rJ, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock );

    //! @brief ... extract a row from the matrix
    //! @param rI ... row
    //! @return row as full matrix
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>  GetRow ( int rI ) const;

    //! @brief ... set a row in this FullMatrix
    //! @param rI ... row
    //! @param rBlock ... new row intries
    void  SetRow ( int rI, const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock );

    //! @brief ... extract a column from this matrix
    //! @param rI ... column
    //! @return column entries as full matrix
    FullVector<T, Eigen::Dynamic>  GetColumn ( int rI ) const;

    //! @brief ... set a column in this matrix
    //! @param rI ... column
    //! @param rBlock ... new column entries
    void  SetColumn ( int rI, const FullVector<T, Eigen::Dynamic>& rBlock );

    //! @brief ... appends columns to this matrix
    //! @param rBlock ... matrix storing the columns
    void  AppendColumns ( const FullVector<T, Eigen::Dynamic>& rBlock );

    //! @brief ... appends columns to this matrix
    //! @param rBlock ... matrix storing the columns
    void  AppendColumns ( const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock );

    //! @brief ... appends rows to this matrix
    //! @param rBlock ... matrix storing the rows
    void  AppendRows ( const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rBlock );

    //! @brief ... calculates the norm of this matrix, i.e. for vectors the Euclidean norm
    //! @return norm of this matrix
    double  Norm() const;

    //! @brief ... computes the sum of all coefficients
    //! @return sum of the entries for each column
    T Sum() const;

    //! @brief ... sets small matrix entries to zero - e.g. for output formating
    //! @param rTolerance ... cut-off tolerance
    void SetSmallEntriesZero(double rTolerance);

    //! @brief ... calculates the sum of the entries for each column
    //! @return sum of the entries for each column
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ColumnwiseSum() const;

    //! @brief ... calculates the sum of the entries for each row
    //! @return sum of the entries for each row
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> RowwiseSum() const;

    //! @brief ... calculates the minimum of the entries for each column
    //! @return the minimum of the entries for each column
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ColumnwiseMinCoeff() const;

    //! @brief ... calculates the minimum of the entries for each row
    //! @return the minimum of the entries for each row
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> RowwiseMinCoeff() const;

    //! @brief ... calculates the maximum of the entries for each column
    //! @return the maximum of the entries for each column
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ColumnwiseMaxCoeff() const;

    //! @brief ... calculates the maximum of the entries for each row
    //! @return the maximum of the entries for each row
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> RowwiseMaxCoeff() const;

    //! @brief ... coefficient wise multiplication of this matrix with another matrix
    //! @param other ... other matrix
    //! @return a matrix which is obtained by a coefficient wise multiplication of this matrix with another matrix
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ElementwiseMul ( const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &other ) const;

    //! @brief ... coefficient wise division of this matrix by another matrix
    //! @param other ... other matrix
    //! @return a matrix which is obtained by a coefficient wise division of this matrix by another matrix
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ElementwiseDiv ( const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &other ) const;

    //! @brief ... coefficient wise reciprocal
    //! @return a matrix which is obtained by coefficient wise reciprocal of this matrix
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> ElementwiseInverse() const;

    //! @brief ... sorts the rows of a matrix based on the entries in column rCol
    FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> SortRow(int rCol);

private:
    //! @brief ... writes a matrix to a file
    //! @param fileStream ... output file stream
    //! @param delimiter_ ... delimiters between the entries in each line
    void WriteToFile ( std::ofstream& fileStream, std::string& delimiter_ ) const;

    //! @brief this is for a workaround to allow partial specialization of single member functions
    template<typename T1, int rows1, int cols1> class params { };

    //! @brief this is the norm member function that is called for anything else than double
    template<typename T1, int rows1, int cols1> double Norm(params<T1, rows1, cols1>)const;

    //! @brief this is the norm member function that is called for doubles
    template<int rows1, int cols1> double Norm(params<double, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for anything else than double
    template<typename T1, int rows1, int cols1> FullMatrix<T,rows,cols> Inverse(params<T1, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for doubles
    template<int rows1, int cols1> FullMatrix<T,rows,cols> Inverse(params<double, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for anything else than double
    template<typename T1, int rows1, int cols1> void EigenValuesSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, params<T1, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for doubles
    template<int rows1, int cols1> void EigenValuesSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, params<double, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for anything else than double
    template<typename T1, int rows1, int cols1> void EigenVectorsSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<T1, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for doubles
    template<int rows1, int cols1> void EigenVectorsSymmetric(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<double, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for anything else than double
    template<typename T1, int rows1, int cols1> void GeneralizedEigenVectorsSymmetric(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rM,FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<T1, rows1, cols1>)const;

    //! @brief this is the inverse member function that is called for doubles
    template<int rows1, int cols1> void GeneralizedEigenVectorsSymmetric(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rM,FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenValues, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rEigenVectors, params<double, rows1, cols1>)const;

    //! @brief this is the SolveCholeskyLapack member function that is called for anything else than double
    template<typename T1, int rows1, int cols1>
    void SolveCholeskyLapack(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRHS, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSolution, params<T1, rows1, cols1>)const;

    //! @brief this is the SolveCholeskyLapack member function that is called for doubles
    template<int rows1, int cols1>
    void SolveCholeskyLapack(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRHS, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSolution, params<double, rows1, cols1>)const;

    //! @brief this is the InverseCholeskyLapack member function that is called for anything else than double
    template<typename T1, int rows1, int cols1>
    void InverseCholeskyLapack(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInverse, params<T1, rows1, cols1>)const;

    //! @brief this is the InverseCholeskyLapack member function that is called for doubles
    template<int rows1, int cols1>
    void InverseCholeskyLapack(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInverse, params<double, rows1, cols1>)const;

}; //class FullMatrix


}//namespace nuto

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>)))

#ifdef ENABLE_SERIALIZATION
//namespace boost
//{
//    //! @brief tell boost how to serialize an Eigen::Matrix
//    template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    inline void serialize(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,const unsigned int file_version)
//    {
//        ar & boost::serialization::make_array(t.data(), t.size());
//    }
//}
#endif // ENABLE_SERIALIZATION
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif // FULL_MATRIX_DEF_H
