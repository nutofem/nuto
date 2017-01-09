// $Id: SparseMatrixCSRVector2Symmetric_Def.h 235 2010-04-22 09:25:38Z arnold2 $

#pragma once

#include "math/SparseMatrixCSRVector2.h"

namespace NuTo
{
template <class T> class SparseMatrixCSRSymmetric;
template <class T> class SparseMatrixCSRVector2General;
//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... class for Symmetric sparse matrices which are stored in CSR format
//!            this file contains the function definitions, whereas the implementation goes either into SparseMatrixCSRVector2Symmetric.h or the cpp file
//!            this was necessary to avoid cyclic inclusions of header templates

template <class T>
class SparseMatrixCSRVector2Symmetric : public SparseMatrixCSRVector2<T>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
    friend class SparseMatrixCSRSymmetric<T>;
    friend class SparseMatrixCSRVector2General<T>;
public:
    //! @brief ... constructor
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns_ ... number of columns
    SparseMatrixCSRVector2Symmetric(int rNumRows_=0, int rNumColumns_=0);

    //! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
    //! @param rFullMatrix ... input matrix (full storage)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
    SparseMatrixCSRVector2Symmetric(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rFullMatrix, double rAbsoluteTolerance = 0, double rRelativeTolerance = 1e-14);

    //! @brief ... create sparse matrix with vector of vector from standard CSR format
    //! @param rCSRMatrix ... input matrix (full storage)
    SparseMatrixCSRVector2Symmetric(const SparseMatrixCSRSymmetric<T>& rCSRMatrix);


    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const override;

    //! @brief ... resize matrix
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns ... number of columns
    void Resize(int rNumRows_, int rNumColumns) override;

    //! @brief ... returns whether the matrix is symmetric or unsymmetric
    //! @return true if the matrix is symmetric and false if the matrix is unsymmetric
    bool IsSymmetric() const override;

    //! @brief ... returns the number of columns
    //! @return number of columns
    int GetNumColumns() const override;

    //! @brief ... add nonzero entry to matrix
    //! @param rRow ... row of the nonzero entry (zero based indexing!!!)
    //! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
    //! @param rValue ... value of the nonzero entry
    void AddValue(int rRow, int rColumn, const T& rValue) override;

    //! @brief ... return the matrix type
    NuTo::eSparseMatrixType GetSparseMatrixType()const override;

    //! @brief ... import matrix from slang object stored in  a text file
    //! @param rFileName ... file name
    void ImportFromSLangText(const char* rFileName) override;

    //! @brief ... write nonzero matrix entries into a matrix
    //! @return ... the matrix
    virtual Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> ConvertToFullMatrix() const override;

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
    friend SparseMatrixCSRVector2Symmetric<T> operator+ (SparseMatrixCSRVector2Symmetric<T> rLhs, const SparseMatrixCSRVector2Symmetric<T> &rRhs )
    {
        rLhs += rRhs;
        return rLhs;
    }

    friend SparseMatrixCSRVector2Symmetric<T> operator- (SparseMatrixCSRVector2Symmetric<T> rLhs, const SparseMatrixCSRVector2Symmetric<T> &rRhs )
    {
        rLhs -= rRhs;
        return rLhs;
    }


    //! @brief ... subtract two matrices
    //! @param rOther ... Symmetric sparse matrix stored in the CSRVector2 format
    //! @return reference to this matrix
    SparseMatrixCSRVector2Symmetric<T>& operator-=  ( const SparseMatrixCSRVector2Symmetric<T> &rOther );

    //! @brief ... add two matrices
    //! @param rOther ... Symmetric sparse matrix stored in the CSRVector2 format
    //! @return reference to this matrix
    SparseMatrixCSRVector2Symmetric<T>& operator+=  ( const SparseMatrixCSRVector2Symmetric<T> &rOther ) override;

    //! @brief ... matrix - matrix multiplication
    //! @param rOther ... Symmetric sparse matrix stored in the CSR format
    //! @return Symmetric sparse matrix stored in the CSR format
    SparseMatrixCSRVector2Symmetric<T> operator* ( const SparseMatrixCSRVector2Symmetric<T> &rOther ) const;

    //! @brief ... multiplies the matrix with an scalar value
    //! @param rOther ... scalar value
    //! @return ... the multiplied matrix (sparse csr storage)
    friend SparseMatrixCSRVector2Symmetric<T> operator* (SparseMatrixCSRVector2Symmetric<T> rLhs,  const T &rRhs )
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
    SparseMatrixCSRVector2Symmetric<T> Transpose() const;

    //! @brief ... calculates the sum of all entries
    T Sum() const override;


    //! @brief ... add the scaled other matrix
    //! @param rOther ... other matrix
    //! @param rFactor ... scalar factor
    void AddScal(const SparseMatrixCSRVector2<T> &rOther, T rFactor) override;


    //! @brief ... add the scaled other matrix
    //! @param rOther ... other matrix
    //! @param rFactor ... scalar factor
    void AddScal(const SparseMatrixCSRSymmetric<T> &rOther, T rFactor);

    //! @brief ... adds \f$\boldsymbol{A}^T\,\boldsymbol{B}\,\boldsymbol{C} \, c\f$ to the matrix
    //! @remark part of \f$\boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KK}\,\boldsymbol{C}_{mat} \, c\f$
    void Add_TransA_B_C_Scal(
            const NuTo::SparseMatrixCSRVector2<T>& rA,
            const NuTo::SparseMatrixCSRVector2<T>& rB,
            const NuTo::SparseMatrixCSRVector2<T>& rC, T rScalar) override;

    //! @brief ... subtract t\f$(\left(\boldsymbol{A}^T\boldsymbol{B} + \boldsymbol{C} \boldsymbol{D}\right)\,c)\f$ from the matrix  with \f$ \boldsymbol{B} = \boldsymbol{C}^T\f$ and \f$ \boldsymbol{D} = \boldsymbol{A}\f$
    //! @remark part of \f$(\boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KJ} + boldsymbol{M}_{JK}\,\boldsymbol{D}_{mat})\,c\f$
    void Sub_TransA_B_Plus_C_D_Scal(
            const SparseMatrixCSRVector2<T>& rA,
            const SparseMatrixCSRVector2<T>& rB,
            const SparseMatrixCSRVector2<T>& rC,
            const SparseMatrixCSRVector2<T>& rD, T rScalar) override;


    //! @brief ... reorder columns of the matrix
    //! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
    void ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering) override;

    //! @brief ... print info about the object
    void Info() const override;

    //! @brief ... returns a random matrix
    //! @param rDimension ... number of rows, number of columns
    //! @param rDensity ... approximate density = numValues / (rNumRows*rNumColumns)
    //! @param rSeed ... random seed
    //! @return random Matrix
    static SparseMatrixCSRVector2Symmetric<T> Random(int rDimension, double rDensity, int rSeed = 0);


    NuTo::SparseMatrixCSRVector2Symmetric<T>& AsSparseMatrixCSRVector2Symmetric()override;
#ifndef SWIG
    const NuTo::SparseMatrixCSRVector2Symmetric<T>& AsSparseMatrixCSRVector2Symmetric()const override;
#endif


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization SparseMatrixCSRVector2Symmetric \n";
#endif
        ar & boost::serialization::make_nvp("SparseMatrixCSRVector2Symmetric",boost::serialization::base_object< SparseMatrixCSRVector2<T> >(*this));
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization SparseMatrixCSRVector2Symmetric \n";
#endif
    }
#endif // ENABLE_SERIALIZATION


protected:
};
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
namespace boost{
template <class T>
struct is_virtual_base_of <NuTo::SparseMatrixCSRVector2<T>, NuTo::SparseMatrixCSRVector2Symmetric<T>>: public mpl:: true_ {};
}
BOOST_CLASS_EXPORT_KEY(NuTo::SparseMatrixCSRVector2Symmetric<double>)
BOOST_CLASS_EXPORT_KEY(NuTo::SparseMatrixCSRVector2Symmetric<int>)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION

