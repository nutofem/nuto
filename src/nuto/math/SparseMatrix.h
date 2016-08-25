// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/Matrix.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{
template <class T, int rows, int cols> class FullMatrix;
template<class T> class SparseMatrixCSRGeneral;
template<class T> class SparseMatrixCSRSymmetric;
template<class T> class SparseMatrixCSRVector2General;
template<class T> class SparseMatrixCSRVector2Symmetric;

enum class eSparseMatrixType;

//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... abstract base class for sparse matrices
template <class T>
class SparseMatrix : public Matrix<T>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief ... constructor
    SparseMatrix() : Matrix<T>()
    {
        this->mOneBasedIndexing=false;
        this->mPositiveDefinite=false;
    }

    SparseMatrix(const SparseMatrix<T>& rOther)
    {
        this->mOneBasedIndexing=rOther.mOneBasedIndexing;
        this->mPositiveDefinite=rOther.mPositiveDefinite;
    }

    virtual ~SparseMatrix() = default;

    //! @brief ... get number of non-zero entries
    //! @return number of non-zero matrix entries
    virtual int GetNumEntries() const = 0;

    //! @brief ... get number of rows
    //! @return number of rows
    virtual int GetNumRows() const = 0;

    //! @brief ... get number of columns
    //! @return number of columns
    virtual int GetNumColumns() const = 0;

    //! @brief ... resize the matrix and initialize everything to zero
    //! @param  rRow ... number of rows
    //! @param  rCol ... number of columns
    virtual void Resize(int rRow, int rCol) = 0;

    //! @brief ... sets all the values to zero while keeping the structure of the matrix constant, this is interesting for stiffness matrices to use the same matrix structure
    virtual void SetZeroEntries() = 0;

	//! @brief ... add nonzero entry to matrix
    //! @param row ... row of the nonzero entry (zero based indexing!!!)
    //! @param column ... column of the nonzero entry (zero based indexing!!!)
    //! @param value ... value of the nonzero entry
    virtual void AddValue(int row, int column, const T& value) = 0;

    //! @brief ... switch to one based indexing
    virtual void SetOneBasedIndexing() = 0;

    //! @brief ... switch to zero based indexing
    virtual void SetZeroBasedIndexing() = 0;

    //! @brief ... get type of indexing
    //! @return true if one based indexing / false if zero based indexing
    inline bool HasOneBasedIndexing() const
    {
        return this->mOneBasedIndexing;
    }

    //! @brief ... get type of indexing
    //! @return false if one based indexing / true if zero based indexing
    inline bool HasZeroBasedIndexing() const
    {
        return ! this->mOneBasedIndexing;
    }

    //! @brief ... get definiteness of matrix
    //! @return true if matrix is positive definite / false if matrix is indefinite
    inline bool IsPositiveDefinite() const
    {
        return this->mPositiveDefinite;
    }

    //! @brief ... set definiteness of matrix to positive definite
    inline void SetPositiveDefinite()
    {
        this->mPositiveDefinite = true;
    }

    //! @brief ... set definiteness of matrix to indefinite
    inline void SetIndefinite()
    {
        this->mPositiveDefinite = false;
    }

    //! @brief ... return the matrix type
    virtual NuTo::eSparseMatrixType GetSparseMatrixType()const=0;

    //! @brief ... symmetry of the matrix
    //! @return ... true if the matrix is symmetric, false otherwise
    virtual bool IsSymmetric() const = 0;

    //! @brief ... remove zero entries from matrix (all entries with an absolute value which is smaller than a prescribed tolerance)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelativeTolerance ... relative tolerance (this value is multiplied with the largest matrix entry (absolute values))
    virtual int RemoveZeroEntries(double rAbsoluteTolerance = 0, double rRelativeTolerance = 0) = 0;

    //! @brief ... returns true if the matrix allows parallel assembly using openmp with maximum independent sets
    //! this is essentially true, if adding a value to a specific row does not change the storage position of values in other rows
    //! until now, this is only true for SparseMatrixCSRVector2
    virtual bool AllowParallelAssemblyUsingMaximumIndependentSets()const
    {
    	return false;
    }

    //! @brief ... multiply sparse matrix with a full matrix
    //! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
    //! @return ... full matrix
    virtual NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> operator* (const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &rMatrix) const = 0;

    //! @brief ... add sparse matrix
    //! @param rMatrix ... sparse matrix
    //! @return ... this
	virtual NuTo::SparseMatrix<T>& operator += (const SparseMatrixCSRSymmetric<T>& rMatrix)
	{
    	throw MathException("[NuTo::SparseMatrix<T>& operator += (const SparseMatrixCSRSymmetric<T> rMatrix)] not implemented for this matrix type.");
	}

	//! @brief ... add sparse matrix
    //! @param rMatrix ... sparse matrix
    //! @return ... this
	virtual NuTo::SparseMatrix<T>& operator += (const SparseMatrixCSRVector2Symmetric<T>& rMatrix)
	{
    	throw MathException("[NuTo::SparseMatrix<T>& operator += (const SparseMatrixCSRVector2Symmetric<T> rMatrix)] not implemented for this matrix type.");
	}

    virtual NuTo::FullMatrix<T,Eigen::Dynamic,Eigen::Dynamic> ConvertToFullMatrixDouble()
    {
        return NuTo::FullMatrix<T,Eigen::Dynamic,Eigen::Dynamic>(*this);
    }

    virtual SparseMatrixCSRGeneral<T>& AsSparseMatrixCSRGeneral()
    {
    	throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRGeneral.");
    }

    virtual SparseMatrixCSRSymmetric<T>& AsSparseMatrixCSRSymmetric()
    {
        throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRSymmetric.");
    }

    virtual SparseMatrixCSRVector2General<T>& AsSparseMatrixCSRVector2General()
    {
        throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRVector2General.");
    }

    virtual SparseMatrixCSRVector2Symmetric<T>& AsSparseMatrixCSRVector2Symmetric()
    {
        throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRVector2Symmetric.");
    }

#ifndef SWIG

    virtual NuTo::FullMatrix<T,Eigen::Dynamic,Eigen::Dynamic> ConvertToFullMatrixDouble() const
    {
        return NuTo::FullMatrix<T,Eigen::Dynamic,Eigen::Dynamic>(*this);
    }

    virtual const SparseMatrixCSRGeneral<T>& AsSparseMatrixCSRGeneral()const
    {
    	throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRGeneral.");
    }

    virtual const SparseMatrixCSRSymmetric<T>& AsSparseMatrixCSRSymmetric()const
    {
        throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRSymmetric.");
    }

    virtual const SparseMatrixCSRVector2General<T>& AsSparseMatrixCSRVector2General()const
    {
        throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRVector2General.");
    }

    virtual const SparseMatrixCSRVector2Symmetric<T>& AsSparseMatrixCSRVector2Symmetric()const
    {
        throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] matrix is not of type SparseMatrixCSRVector2Symmetric.");
    }
#endif // SWIG






#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of SparseMatrix" << std::endl;
#endif
        ar & boost::serialization::make_nvp("Matrix",boost::serialization::base_object< Matrix<T> >(*this));
        ar & BOOST_SERIALIZATION_NVP(mOneBasedIndexing)
           & BOOST_SERIALIZATION_NVP(mPositiveDefinite);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of SparseMatrix" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

    //! @brief ... write non-zero matrix entries into a matrix
    //! @param rMatrix ... the matrix
    virtual void WriteEntriesToMatrix(NuTo::Matrix<T>& rMatrix) const = 0;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId() const;

    //! @brief Calculate the largest matrix entry
    //! @param rResultOutput ... largest matrix entry
    T Max() const
    {
        T max;
        int row, col;
        MaxEntry(row, col, max);
        return max;
    }

    //! @brief Calculate the smallest matrix entry
    //! @param rResultOutput ... smallest matrix entry
    T Min() const
    {
        T min;
        int row, col;
        MinEntry(row, col, min);
        return min;
    }

    virtual void MaxEntry(int& row_output, int& column_output, T& result_output) const = 0;

    virtual void MinEntry(int& row_output, int& column_output, T& result_output) const = 0;

    //! @brief calculates the largest absolute matrix entry with the corresponding position
    //! @param rRow ... row
    //! @param rCol ... column
    //! @return ... largest absolute matrix entry
    T AbsMax(int& rRow, int& rCol) const
    {
        T max, min;
        int row, col;

        MaxEntry(rRow, rCol, max);
        MinEntry( row,  col, min);

        if (std::abs(max) > std::abs(min))
        {
            // (max, rRow, rCol) is abs maximum
            return max;
        }
        else
        {
            // (min, row, col) is abs maximum
            rRow = row;
            rCol = col;
            return min;
        }
    }

    //! @brief calculates the largest absolute matrix entry
    //! @return ... largest absolute matrix entry
    T AbsMax() const
    {
        int row, col;
        return AbsMax(row, col);
    }

protected:
    //! @brief ... internal indexing of the matrix (true if one based indexing / false if zero based indexing)
    bool mOneBasedIndexing;
    //! @brief ... definiteness of the matrix (true if positive definite / false if indefinite)
    bool mPositiveDefinite;
};


}
