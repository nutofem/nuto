// $Id: SparseMatrixCSRVector2.h 235 2010-04-22 09:25:38Z arnold2 $

#ifndef SPARSE_MATRIX_CSR_VECTOR2_H
#define SPARSE_MATRIX_CSR_VECTOR2_H
#include <vector>
#include <string>
#include <climits>
#include <assert.h>
#include <iostream>
#include <fstream>  //for file acces

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION
#include <boost/foreach.hpp>

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
template <class T> class SparseMatrixCSRGeneral;
//! @author Jörg F. Unger, NU
//! @date July 2010
//! @brief ... abstract base class for sparse matrices which are stored in CSR format
//! with a vector in vector format in order to reduce the amount required for inserting to that matrix when building the stiffness matrix
template <class T>
class SparseMatrixCSRVector2 : public SparseMatrix<T>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
    friend class NuTo::SparseMatrixCSRGeneral<T>;

public:
    //! @brief ... constructor
    //! @param rNumRows_ ... number of rows
    //! @param rNumReserveEntries_ ... number of entries per row for which memory is reserved (optional)
    SparseMatrixCSRVector2(int rNumRows_) : SparseMatrix<T>()
    {
        // check for overflow
        assert(rNumRows_ < INT_MAX);
        assert(rNumRows_ >= 0);
        this->mValues.resize(rNumRows_);
        this->mColumns.resize(rNumRows_);
    }

    //! @brief ... reserve memory for non-zero matrix entries
    //! @param rNumReserveEntries_ ... number of entries per row for which memory is reserved
    void Reserve(unsigned int rNumReserveEntries_)
    {
        for (unsigned int count=0; count<this->mValues.size(); count++)
        {
        	this->mValues[count].reserve(rNumReserveEntries_);
        	this->mColumns[count].reserve(rNumReserveEntries_);
        }
    }

    //! @brief ... returns the number of non-zero matrix entries
    //! @return number of non-zero matrix entries
    int GetNumEntries() const
    {
        unsigned int counter(0);
        for (unsigned int count=0; count<this->mValues.size(); count++)
        {
         	counter += this->mValues[count].size();
        }
    	return counter;
    }

    //! @brief ... returns the number of rows
    //! @return number of rows
    int GetNumRows() const
    {
        return this->mValues.size();
    }

    //! @brief ... switch to one based indexing (only internal indexing, interface still uses zero based indexing)
    void SetOneBasedIndexing()
    {
        if (! this->mOneBasedIndexing)
        {
            for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
            {
                for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
                {
                	this->mColumns[row_count][col_count]++;
                }
            }
            this->mOneBasedIndexing = true;
        }
    }

    //! @brief ... switch to zero based indexing (only internal indexing, interface still uses zero based indexing)
    void SetZeroBasedIndexing()
    {
        if (this->mOneBasedIndexing)
        {
            for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
            {
                for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
                {
                	this->mColumns[row_count][col_count]--;
                }
            }
            this->mOneBasedIndexing = false;
        }
    }


    //! @brief ... write non-zero matrix entries into a matrix
    //! @param rFullMatrix ... the full matrix
    virtual void WriteEntriesToMatrix(NuTo::Matrix<T>& rMatrix) const = 0;

    //! @brief ... import matrix from slang object stored in  a text file
    //! @param rFileName ... file name
    virtual void ImportFromSLangText(const char* rFileName) = 0;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId() const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp("SparseMatrix",boost::serialization::base_object< SparseMatrix<T> >(*this));
        ar & BOOST_SERIALIZATION_NVP(mColumns)
           & BOOST_SERIALIZATION_NVP(mValues);
   }
#endif  // ENABLE_SERIALIZATION

    //! @brief performs a monadic operator on all matrix entries
    //! @param rMOperator        Monadic Operator
    void Map(const NuTo::MonadicOperator<T>* rMOperator)
    {
        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
                mValues[row_count][col_count] = rMOperator->Evaluate(mValues[row_count][col_count]);
            }
        }
    }

    //! @brief performs a dyadic operator on all matrix entries with another given value
    //! @param rDOperator        Dyadic Operator
    //! @param rValue ... value
    void Map(const NuTo::DyadicOperator<T>* rDOperator, const T& rValue)
    {
        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
                mValues[row_count][col_count] = rDOperator->Evaluate(mValues[row_count][col_count],rValue);
            }
        }
    }

    //! @brief ... multiplies the matrix with an scalar value
    //! @param rOther ... scalar value
    //! @return ... the multiplied matrix (sparse csr storage)
    SparseMatrixCSRVector2<T>& operator*=  ( const T &rOther )
    {
        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
                mValues[row_count][col_count] *= rOther;
            }
        }
        return *this;
    }

    //! @brief Calculate the largest matrix entry
    //! @param rResultOutput ... largest matrix entry
    virtual T Max()
    {
        T max;
        Max(max);
        return max;
    }


#ifndef SWIG
    //! @brief Determine the largest matrix entry
    //! @param rResultOutput ... largest matrix entry
    void Max(T& rResultOutput)
    {
        //search for first entry
        bool first_entry_found(false);
        for (unsigned int row_count = 0; row_count < this->mColumns.size() && first_entry_found==false; row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
           		rResultOutput = mValues[row_count][col_count];
           		first_entry_found=true;
            }
        }

        if (first_entry_found==false)
            throw MathException("[NuTo::SparseMatrixCSRVector2::Max] Maximum for matrix with zero entries cannot be calculated.");

        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
            	if (mValues[row_count][col_count]>rResultOutput)
            		rResultOutput = mValues[row_count][col_count];
            }
        }
    }
#endif

    //! @brief Calculate the smallest matrix entry
    //! @param rResultOutput ... smallest matrix entry
    virtual T Min()
    {
        T min;
        Min(min);
        return min;
    }

#ifndef SWIG
    //! @brief Calculate the smallest matrix entry
    //! @param rResultOutput ... smallest matrix entry
    virtual void Min(T& rResultOutput)
    {
        //search for first entry
        bool first_entry_found(false);
        for (unsigned int row_count = 0; row_count < this->mColumns.size() && first_entry_found==false; row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
           		rResultOutput = mValues[row_count][col_count];
           		first_entry_found=true;
            }
        }

        if (first_entry_found==false)
            throw MathException("[NuTo::SparseMatrixCSRVector2::Min] Minimum for matrix with zero entries cannot be calculated.");

        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
            	if (mValues[row_count][col_count]<rResultOutput)
            		rResultOutput = mValues[row_count][col_count];
            }
        }
    }
#endif

    //! @brief Calculate the maximum matrix entry
    //! @param rRowOutput ... row index of the maximum matrix entry
    //! @param rColumnOutput ... column index of the maximum matrix entry
    //! @param rResultOutput ... maximum matrix entry
    virtual void Max(int& rRowOutput, int& rColumnOutput, T& rResultOutput)
    {
        //search for first entry
        bool first_entry_found(false);
        for (unsigned int row_count = 0; row_count < this->mColumns.size() && first_entry_found==false; row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
           		rResultOutput = mValues[row_count][col_count];
           		if (this->mOneBasedIndexing)
           		    rRowOutput = row_count+1;
           		else
           		    rRowOutput = row_count;

       		    rColumnOutput = mColumns[row_count][col_count];
          		first_entry_found=true;
            }
        }

        if (first_entry_found==false)
            throw MathException("[NuTo::SparseMatrixCSRVector2::Max] Maximum for matrix with zero entries cannot be calculated.");

        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
            	if (mValues[row_count][col_count]>rResultOutput)
            	{
            		rResultOutput = mValues[row_count][col_count];
					if (this->mOneBasedIndexing)
						rRowOutput = row_count+1;
					else
						rRowOutput = row_count;

					rColumnOutput = mColumns[row_count][col_count];
            	}
            }
        }
    }

    //! @brief Calculate the minimum matrix entry
    //! @param rRowOutput ... row index of the minimum matrix entry
    //! @param rColumnOutput ... column index of the minimum matrix entry
    //! @param rResultOutput ... minimum matrix entry
    virtual void Min(int& rRowOutput, int& rColumnOutput, T& rResultOutput)
    {
        //search for first entry
        bool first_entry_found(false);
        for (unsigned int row_count = 0; row_count < this->mColumns.size() && first_entry_found==false; row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
           		rResultOutput = mValues[row_count][col_count];
           		if (this->mOneBasedIndexing)
           		    rRowOutput = row_count+1;
           		else
           		    rRowOutput = row_count;

       		    rColumnOutput = mColumns[row_count][col_count];
          		first_entry_found=true;
            }
        }

        if (first_entry_found==false)
            throw MathException("[NuTo::SparseMatrixCSRVector2::Min] Minimum for matrix with zero entries cannot be calculated.");

        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
            	if (mValues[row_count][col_count]<rResultOutput)
            	{
            		rResultOutput = mValues[row_count][col_count];
					if (this->mOneBasedIndexing)
						rRowOutput = row_count+1;
					else
						rRowOutput = row_count;

					rColumnOutput = mColumns[row_count][col_count];
            	}
            }
        }
    }

    //! @brief ... sets all the values to zero while keeping the structure of the matrix constant, this is interesting for stiffness matrices to use the same matrix structure
    void SetZeroEntries()
    {
        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
            	mValues[row_count][col_count]= 0;
            }
        }
    }

    //! @brief ... remove zero entries from matrix (all entries with an absolute value which is smaller than a prescribed tolerance)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelativeTolerance ... relative tolerance (this value is multiplied with the largest matrix entry (absolute values))
    int RemoveZeroEntries(double rAbsoluteTolerance = 0, double rRelativeTolerance = 0);

    // @brief ... returns the data values to allow other methods to loop over the nonzero entries
    // @return mValues
    const std::vector<std::vector<T> >& GetValues()const
    {
        return  mValues;
    }

    // @brief ... returns the data values to allow other methods to loop over the nonzero entries
    // @return mValues
    const std::vector<std::vector<int> >& GetColumns()const
    {
        return  mColumns;
    }

    //! @brief ... returns true if the matrix allows parallel assembly using openmp with maximum independent sets
    //! this is essentially true, if adding a value to a specific row does not change the storage position of values in other rows
    //! until now, this is only true for SparseMatrixCSRVector2
    bool AllowParallelAssemblyUsingMaximumIndependentSets()const
    {
    	return true;
    }

protected:
    //! @brief value of nonzero matrix entries
    std::vector<std::vector<T> > mValues;
    //! @brief columns of nonzero matrix entries
    std::vector<std::vector<int> > mColumns;
};
}
#endif // SPARSE_MATRIX_CSR_VECTOR2_H
