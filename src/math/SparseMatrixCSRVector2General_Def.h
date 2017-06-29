#pragma once

#include "base/Exception.h"
#include "math/SparseMatrixCSRVector2.h"

namespace NuTo
{
template <class T> class SparseMatrixCSRGeneral;
template <class T> class SparseMatrixCSRVector2Symmetric;

//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... class for general sparse matrices which are stored in CSR format
//!            this file contains the function definitions, whereas the implementation goes either into SparseMatrixCSRVector2General.h or the cpp file
//!            this was necessary to avoid cyclic inclusions of header templates

template <class T>
class SparseMatrixCSRVector2General : public SparseMatrixCSRVector2<T>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
    friend class SparseMatrixCSRVector2Symmetric<T>;
public:
    //! @brief ... constructor
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns_ ... number of columns
    SparseMatrixCSRVector2General(int rNumRows=0, int rNumColumns=0);

    //! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
    //! @param rFullMatrix ... input matrix (full storage)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
    SparseMatrixCSRVector2General(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rFullMatrix, double rAbsoluteTolerance = 0, double rRelativeTolerance = 1e-14);

    //! @brief ... create sparse matrix with vector of vector from standard CSR format
    //! @param rCSRMatrix ... input matrix (full storage)
    SparseMatrixCSRVector2General(const SparseMatrixCSRGeneral<T>& rCSRMatrix);

    //! @brief ... create sparse matrix with vector of vector from symmetric sparse matrix
    //! @param rCSRVector2Symmetric ... input matrix (full storage)
    SparseMatrixCSRVector2General(const SparseMatrixCSRVector2Symmetric<T>& rCSRVector2Symmetric);

    //! @brief ... resize matrix
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns ... number of columns
    void Resize(int rNumRows, int rNumColumns) override;

    //! @brief ... returns whether the matrix is symmetric or unsymmetric
    //! @return true if the matrix is symmetric and false if the matrix is unsymmetric
    bool IsSymmetric() const override;

    //! @brief ... returns the number of columns
    //! @return number of columns
    int GetNumColumns() const override;

    void Info() const override;

    //! @brief ... add nonzero entry to matrix
    //! @param rRow ... row of the nonzero entry (zero based indexing!!!)
    //! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
    //! @param rValue ... value of the nonzero entry
    void AddValue(int rRow, int rColumn, const T& rValue) override;

    //! @brief ... return the matrix type
    NuTo::eSparseMatrixType GetSparseMatrixType()const override;

    //! @brief ... write nonzero matrix entries into a matrix
    //! @return ... the matrix
    virtual Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> ConvertToFullMatrix() const override;

    //! @brief ... returns the symmetric part of the matrix 0.5*(A+A^T)
    //! @return symmetric part
    SparseMatrixCSRVector2Symmetric<T> SymmetricPart() const;

#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const;

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType);
#endif // ENABLE_SERIALIZATION

#ifndef SWIG
    friend SparseMatrixCSRVector2General<T> operator+ (SparseMatrixCSRVector2General<T> rLhs, const SparseMatrixCSRVector2General<T> &rRhs )
    {
        rLhs += rRhs;
        return rLhs;
    }

    friend SparseMatrixCSRVector2General<T> operator+ (SparseMatrixCSRVector2General<T> rLhs, const SparseMatrixCSRVector2Symmetric<T> &rRhs )
    {
        rLhs += rRhs;
        return rLhs;
    }

    friend SparseMatrixCSRVector2General<T> operator- (SparseMatrixCSRVector2General<T> rLhs, const SparseMatrixCSRVector2General<T> &rRhs )
    {
        rLhs -= rRhs;
        return rLhs;
    }

    friend SparseMatrixCSRVector2General<T> operator- (SparseMatrixCSRVector2General<T> rLhs, const SparseMatrixCSRVector2Symmetric<T> &rRhs )
    {
        rLhs -= rRhs;
        return rLhs;
    }

    //! @brief ... subtract two matrices
    //! @param rOther ... general sparse matrix stored in the CSRVector2 format
    //! @return general sparse matrix stored in the CSR format
//    SparseMatrixCSRVector2General<T> operator- ( const SparseMatrixCSRVector2General<T> &rOther );

    //! @brief ... subtract two matrices
    //! @param rOther ... general sparse matrix stored in the CSRVector2 format
    //! @return reference to this matrix
    SparseMatrixCSRVector2General<T>& operator-=  ( const SparseMatrixCSRVector2General<T> &rOther );

    //! @brief ... add two matrices
    //! @param rOther ... general sparse matrix stored in the CSRVector2 format
    //! @return reference to this matrix
    SparseMatrixCSRVector2General<T>& operator+=  ( const SparseMatrixCSRVector2General<T> &rOther );

    //! @brief ... add two matrices
    //! @param rOther ... symmetric sparse matrix stored in the CSRVector2 format
    //! @return reference to this matrix
    SparseMatrixCSRVector2General<T>& operator+=  ( const SparseMatrixCSRVector2Symmetric<T> &rOther ) override;

    //! @brief ... matrix - matrix multiplication
    //! @param rOther ... general sparse matrix stored in the CSR format
    //! @return general sparse matrix stored in the CSR format
    SparseMatrixCSRVector2General<T> operator* ( const SparseMatrixCSRVector2General<T> &rOther ) const;

    //! @brief ... multiplies the matrix with an scalar value
    //! @param rOther ... scalar value
    //! @return ... the multiplied matrix (sparse csr storage)
    friend SparseMatrixCSRVector2General<T> operator* (SparseMatrixCSRVector2General<T> rLhs,  const T &rRhs )
    {
        rLhs *= rRhs;
        return rLhs;
    }

#endif // SWIG

    //! @brief ... multiply sparse matrix with a full matrix
    //! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
    //! @return ... full matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> operator* (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &rMatrix) const override;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TransMult(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const override;

    //! @brief ... calculate the transpose of the matrix (transpose row and columns)
    //! @return ... transpose of this matrix (sparse csr storage)
    SparseMatrixCSRVector2General<T> Transpose() const;

    //! @brief ... calculates the sum of all entries
    T Sum() const override;

    //! @brief ... add the scaled other matrix
    //! @param rOther ... other matrix
    //! @param rFactor ... scalar factor
    void AddScal(const SparseMatrixCSRVector2<T> &rOther, T rFactor) override;

    //! @brief ... add the scaled other matrix
    //! @param rOther ... other matrix
    //! @param rFactor ... scalar factor
    void AddScal(const SparseMatrixCSRVector2Symmetric<T> &rOther, T rFactor);

    //! @brief ... add the scaled other matrix
    //! @param rOther ... other matrix
    //! @param rFactor ... scalar factor
    void AddScal(const SparseMatrixCSRGeneral<T> &rOther, T rFactor);

    //! @brief ... add the scaled other matrix
    //! @param rOther ... other matrix
    //! @param rFactor ... scalar factor
    void AddScal(const SparseMatrixCSRSymmetric<T> &rOther, T rFactor);

    //! @brief ... adds \f$\boldsymbol{A}^T\,\boldsymbol{B}\,\boldsymbol{C} \, c\f$ to the matrix
    //! @remark part of \f$\boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KK}\,\boldsymbol{C}_{mat} \, c\f$
    void Add_TransA_B_C_Scal(
            const SparseMatrixCSRVector2<T>& rA,
            const SparseMatrixCSRVector2<T>& rB,
            const SparseMatrixCSRVector2<T>& rC, T rScalar) override;

    //! @brief ... subtract t\f$(\left(\boldsymbol{A}^T\boldsymbol{B}^T + \boldsymbol{B} \boldsymbol{D}\right)\,c)\f$ from the matrix
    //! @remark part of \f$(\boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KJ} + \boldsymbol{M}_{JK}\,\boldsymbol{C}_{mat})\,c\f$
    void Sub_TransA_B_Plus_C_D_Scal(
            const SparseMatrixCSRVector2<T>& rA,
            const SparseMatrixCSRVector2<T>& rB,
            const SparseMatrixCSRVector2<T>& rC,
            const SparseMatrixCSRVector2<T>& rD, T rScalar) override;

    //! @brief ... perform Gauss algorithm (matrix and right hand side are reordered and modified)
    //! @param rRhs ... right-hand side vector (input and output object)
    //! @param rMappingNewToInitialOrdering ... mapping from new ordering to initial ordering (output object)
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to new ordering (output object)
    //! @param rRelativeTolerance ... relative tolerance for zero matrix entries
    void Gauss(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance = 1e-14);

    //! @brief ... perform Gauss algorithm (matrix and right hand side are reordered and modified)
    //! @param rRhs ... right-hand side vector (input and output object)
    //! @param rMappingNewToInitialOrdering ... mapping from new ordering to initial ordering (output object)
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to new ordering (output object)
    //! @param rRelativeTolerance ... relative tolerance for zero matrix entries
    //! @remark converts everything from CSRVector2 to CSR, calculates Gauss, converts everything back to CSRVector2. Fix that.
    void Gauss(NuTo::SparseMatrixCSRVector2<T>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance = 1e-14) override;

    //! @brief ... reorder columns of the matrix
    //! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
    void ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering) override;

    //! @brief ... remove columns from the end of the matrix
    //! @param rNumColumn ... number of colums to be removed
    void RemoveLastColumns(unsigned int rNumColumns) override;

    //! @brief ... Concatenate columns from another matrix
    //! @param rOther ... other matrix with same number of rows
    void ConcatenateColumns(const SparseMatrixCSRVector2General<T>& rOther);

    //! @brief ... Concatenate rows from another matrix
    //! @param rOther ... other matrix with same number of columns
    void ConcatenateRows(const SparseMatrixCSRVector2General<T>& rOther);

    //! @brief ... returns a random matrix
    //! @param rNumRows ... number of rows
    //! @param rNumColumns ... number of columns
    //! @param rDensity ... approximate density = numValues / (rNumRows*rNumColumns)
    //! @param rSeed ... random seed
    //! @return random Matrix
    static SparseMatrixCSRVector2General<T> Random(int rNumRows, int rNumColumns, double rDensity, int rSeed = 0);

    NuTo::SparseMatrixCSRVector2General<T>& AsSparseMatrixCSRVector2General()override;
#ifndef SWIG
    const NuTo::SparseMatrixCSRVector2General<T>& AsSparseMatrixCSRVector2General()const override;
#endif

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization SparseMatrixCSRVector2General \n";
#endif
        ar & boost::serialization::make_nvp("SparseMatrixCSRVector2General",boost::serialization::base_object< SparseMatrixCSRVector2<T> >(*this));
        ar & BOOST_SERIALIZATION_NVP(mNumColumns);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization SparseMatrixCSRVector2General \n";
#endif
    }
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief ... number of columns
    int mNumColumns;

};
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
namespace boost{
template <class T>
struct is_virtual_base_of <NuTo::SparseMatrixCSRVector2<T>, NuTo::SparseMatrixCSRVector2General<T>>: public mpl:: true_ {};
}
BOOST_CLASS_EXPORT_KEY(NuTo::SparseMatrixCSRVector2General<double>)
BOOST_CLASS_EXPORT_KEY(NuTo::SparseMatrixCSRVector2General<int>)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION

