#include <string>
#include "math/SparseMatrix.h"
#include "math/SparseMatrixCSR.h"
#include "math/MathException.h"

namespace NuTo
{

template<>
int SparseMatrixCSR<int>::RemoveZeroEntries(double rAbsoluteTolerance, double rRelativeTolerance)
{
	   double tolerance =  rAbsoluteTolerance;
	    if (rRelativeTolerance > 0)
	    {
	        double maxValue = 0;
	        for (unsigned int entryCount = 0; entryCount < this->mValues.size(); entryCount++)
	        {
	            if (std::abs(this->mValues[entryCount]) > maxValue)
	            {
	                maxValue = std::abs(this->mValues[entryCount]);
	            }
	        }
	        tolerance += rRelativeTolerance * maxValue;
	    }

	    int numRemoved(0);
	    int newPos = 0;
	    int start = newPos;
	    if (this->mOneBasedIndexing)
	    {
	        for (int row = 0; row < this->GetNumRows(); row++)
	        {
	            for (int pos = start; pos < this->mRowIndex[row + 1] -1 ; pos++)
	            {
	                if (std::abs(this->mValues[pos]) > tolerance)
	                {
	                    this->mValues[newPos] = this->mValues[pos];
	                    this->mColumns[newPos] = this->mColumns[pos];
	                    newPos++;
	                }
	                else
	                	numRemoved++;
	            }
	            start = this->mRowIndex[row + 1] - 1;
	            this->mRowIndex[row + 1] = newPos + 1;

	        }
	    }
	    else
	    {
	        for (int row = 0; row < this->GetNumRows(); row++)
	        {
	            for (int pos = start; pos < this->mRowIndex[row + 1]; pos++)
	            {
	                if (std::abs(this->mValues[pos]) > tolerance)
	                {
	                    this->mValues[newPos] = this->mValues[pos];
	                    this->mColumns[newPos] = this->mColumns[pos];
	                    newPos++;
	                }
	                else
	                 	numRemoved++;
	            }
	            start = this->mRowIndex[row + 1];
	            this->mRowIndex[row + 1] = newPos;
	        }
	    }
	    mValues.resize(mValues.size()-numRemoved);
	    mColumns.resize(mColumns.size()-numRemoved);

	    return numRemoved;
}
template<>
int SparseMatrixCSR<double>::RemoveZeroEntries(double rAbsoluteTolerance, double rRelativeTolerance)
{
    double tolerance =  rAbsoluteTolerance;
    if (rRelativeTolerance > 0)
    {
        double maxValue = 0;
        for (unsigned int entryCount = 0; entryCount < this->mValues.size(); entryCount++)
        {
            if (std::abs(this->mValues[entryCount]) > maxValue)
            {
                maxValue = std::abs(this->mValues[entryCount]);
            }
        }
        tolerance += rRelativeTolerance * maxValue;
    }

    int numRemoved(0);
    int newPos = 0;
    int start = newPos;
    if (this->mOneBasedIndexing)
    {
        for (int row = 0; row < this->GetNumRows(); row++)
        {
            for (int pos = start; pos < this->mRowIndex[row + 1] -1 ; pos++)
            {
                if (std::abs(this->mValues[pos]) > tolerance)
                {
                    this->mValues[newPos] = this->mValues[pos];
                    this->mColumns[newPos] = this->mColumns[pos];
                    newPos++;
                }
                else
                	numRemoved++;
            }
            start = this->mRowIndex[row + 1] - 1;
            this->mRowIndex[row + 1] = newPos + 1;

        }
    }
    else
    {
        for (int row = 0; row < this->GetNumRows(); row++)
        {
            for (int pos = start; pos < this->mRowIndex[row + 1]; pos++)
            {
                if (std::abs(this->mValues[pos]) > tolerance)
                {
                    this->mValues[newPos] = this->mValues[pos];
                    this->mColumns[newPos] = this->mColumns[pos];
                    newPos++;
                }
                else
                 	numRemoved++;
            }
            start = this->mRowIndex[row + 1];
            this->mRowIndex[row + 1] = newPos;
        }
    }
    mValues.resize(mValues.size()-numRemoved);
    mColumns.resize(mColumns.size()-numRemoved);

    return numRemoved;
}

}
