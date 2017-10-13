// $Id: SparseMatrixCSRVector2.h 235 2010-04-22 09:25:38Z arnold2 $

#pragma once
#include <vector>
#include <string>
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

#include "nuto/math/FullVector_Def.h"
#include "nuto/math/SparseMatrix.h"

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

    virtual ~SparseMatrixCSRVector2() {}

    SparseMatrixCSRVector2(const SparseMatrixCSRVector2&) = default;

#ifndef SWIG

    SparseMatrixCSRVector2(SparseMatrixCSRVector2&&) = default;

    //! @brief ... copy assignment operator that (in contrast to the default implementation) resizes (*this) first
    //! @param rOther ... other matrix
    SparseMatrixCSRVector2& operator= (const SparseMatrixCSRVector2& rOther)
    {
        this->Resize(rOther.GetNumRows(), rOther.GetNumColumns());
        mValues = rOther.mValues;
        mColumns = rOther.mColumns;
        return *this;
    }


    //! @brief ... move assignment operator that (in contrast to the default implementation) resizes (*this) first
    //! @param rOther ... other matrix
    SparseMatrixCSRVector2& operator= (SparseMatrixCSRVector2&& rOther)
    {
        this->Resize(rOther.GetNumRows(), rOther.GetNumColumns());
        mValues = std::move(rOther.mValues);
        mColumns = std::move(rOther.mColumns);
        return *this;
    }
#endif // SWIG

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
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of SparseMatrixCSRVector2" << std::endl;
#endif
        ar & boost::serialization::make_nvp("SparseMatrixCSRVector2",boost::serialization::base_object< SparseMatrix<T> >(*this));
        ar & BOOST_SERIALIZATION_NVP(mColumns)
           & BOOST_SERIALIZATION_NVP(mValues);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of SparseMatrixCSRVector2" << std::endl;
#endif
   }
#endif // SWIG
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

    //! @brief Calculate the maximum matrix entry
    //! @param rRowOutput ... row index of the maximum matrix entry
    //! @param rColumnOutput ... column index of the maximum matrix entry
    //! @param rResultOutput ... maximum matrix entry
    virtual void MaxEntry(int& rRowOutput, int& rColumnOutput, T& rResultOutput) const override
    {
        rResultOutput = GetFirstValue(rRowOutput, rColumnOutput);

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
    virtual void MinEntry(int& rRowOutput, int& rColumnOutput, T& rResultOutput) const override
    {
        rResultOutput = GetFirstValue(rRowOutput, rColumnOutput);

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
    int RemoveZeroEntries(double rAbsoluteTolerance = 0, double rRelativeTolerance = 0)
    {
        double tolerance =  rAbsoluteTolerance;
         if (rRelativeTolerance > 0)
         {
             double maxValue = this->Max();
             double minValue = this->Min();

             if (std::abs(maxValue)>std::abs(minValue))
                 tolerance += rRelativeTolerance * maxValue;
         }

         int numRemoved(0);
         for (int row = 0; row < this->GetNumRows(); row++)
         {
             unsigned int numRemovedPerRow(0);
             unsigned int newPos(0);
             std::vector<T>& thisValueVec = this->mValues[row];
             std::vector<int>& thisColumnVec = this->mColumns[row];
             for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
             {
                 if (std::abs(thisValueVec[pos]) > tolerance)
                 {
                     thisValueVec[newPos] = thisValueVec[pos];
                     thisColumnVec[newPos] = thisColumnVec[pos];
                     newPos++;
                 }
                 else
                 {
                     numRemovedPerRow++;
                 }
             }
             thisValueVec.resize(thisValueVec.size()-numRemovedPerRow);
             thisColumnVec.resize(thisColumnVec.size()-numRemovedPerRow);
             numRemoved+=numRemovedPerRow;
         }

         return numRemoved;
    }

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

    // @brief ... returns the data values to allow other methods to loop over the nonzero entries
    // @return mValues
    std::vector<std::vector<T> >& GetValuesReference()
    {
        return  mValues;
    }

    // @brief ... returns the data values to allow other methods to loop over the nonzero entries
    // @return mValues
    std::vector<std::vector<int> >& GetColumnsReference()
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


    virtual void AddScal(const SparseMatrixCSRVector2<T>& rOther, T rScalar) = 0;

    //! @brief ... calculates the sum of all entries
    virtual T Sum() const = 0;

    //! @brief ... add the scaled other diagonal matrix
    //! @param rOther ... other vector interpreted as diagonal matrix
    //! @param rFactor ... scalar factor
    void AddScalDiag(const NuTo::FullVector<T, Eigen::Dynamic>& rOther, T rScalar)
    {
        if ((this->GetNumColumns() != rOther.GetNumRows()) || (this->GetNumRows() != rOther.GetNumRows()))
            throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");

        for (int row = 0; row < rOther.GetNumRows(); row++)
            this->AddValue(row, row, rScalar*rOther[row]);
    }


    //! @brief ... adds \f$\boldsymbol{A}^T\,\boldsymbol{B}\,\boldsymbol{C} \, c\f$ to the matrix
    //! @remark part of \f$\boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KK}\,\boldsymbol{C}_{mat} \, c\f$
    //! @param rA ... only Vector2General matrices are reasonable
    //! @param rB ... result is Vector2Symmetric of rB is Vector2Symmetric, result is Vector2General if rB is Vector2General
    //! @param rC ... only Vector2General matrices are reasonable
    virtual void Add_TransA_B_C_Scal(
            const SparseMatrixCSRVector2& rA,
            const SparseMatrixCSRVector2& rB,
            const SparseMatrixCSRVector2& rC, T rScalar) = 0;

    //! @brief ... subtract t\f$(\left(\boldsymbol{A}^T\boldsymbol{B}^T + \boldsymbol{C} \boldsymbol{D}\right)\,c)\f$ from the matrix
    //! @remark part of \f$(\boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KJ} + boldsymbol{M}_{JK}\,\boldsymbol{C}_{mat})\,c\f$
    virtual void Sub_TransA_B_Plus_C_D_Scal(
            const SparseMatrixCSRVector2& rA,
            const SparseMatrixCSRVector2& rB,
            const SparseMatrixCSRVector2& rC,
            const SparseMatrixCSRVector2& rD, T rScalar) = 0;

    //! @brief ... calculates this.Transpose
    virtual NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> TransMult(const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const
    {
        throw NuTo::MathException(std::string("[") + __PRETTY_FUNCTION__ + "] not implemented.");
    }

    //! @brief inverts the matrix coefficient-wise
    void CwiseInvert()
    {
        for (auto& row : mValues)
            for (T& value : row)
                value = 1./value;
    }


    //! @brief ... perform Gauss algorithm (matrix and right hand side are reordered and modified)
    //! @param rRhs ... right-hand side vector (input and output object)
    //! @param rMappingNewToInitialOrdering ... mapping from new ordering to initial ordering (output object)
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to new ordering (output object)
    //! @param rRelativeTolerance ... relative tolerance for zero matrix entries
    virtual void Gauss(NuTo::SparseMatrixCSRVector2<T>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance = 1e-14)
    {
        throw NuTo::MathException(std::string("[") + __PRETTY_FUNCTION__ + "] not implemented.");
    }

    //! @brief ... reorder columns of the matrix
    //! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
    virtual void ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering)
    {
        throw NuTo::MathException(std::string("[") + __PRETTY_FUNCTION__ + "] not implemented.");
    }

    //! @brief ... remove columns from the end of the matrix
    //! @param rNumColumn ... number of colums to be removed
    virtual void RemoveLastColumns(unsigned int rNumColumns)
    {
        throw NuTo::MathException(std::string("[") + __PRETTY_FUNCTION__ + "] not implemented.");
    }

protected:
    //! @brief value of nonzero matrix entries
    std::vector<std::vector<T> > mValues;
    //! @brief columns of nonzero matrix entries
    std::vector<std::vector<int> > mColumns;


private:
    //! @brief returns any value in the matrix
    T GetFirstValue(int& rRow, int& rCol) const
    {
        for (unsigned int row_count = 0; row_count < this->mColumns.size(); row_count++)
        {
            for (unsigned int col_count = 0; col_count < this->mColumns[row_count].size(); col_count++)
            {
                rCol = mColumns[row_count][col_count];
                if (this->mOneBasedIndexing)
                    rRow = row_count+1;
                else
                    rRow = row_count;
                return mValues[row_count][col_count];
            }
        }
        throw MathException("[GetFirstValue] Matrix has no entries.");
    }

};
}
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::SparseMatrixCSRVector2<double>)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::SparseMatrixCSRVector2<int>)
#endif  // SWIG
#endif  // ENABLE_SERIALIZATION
