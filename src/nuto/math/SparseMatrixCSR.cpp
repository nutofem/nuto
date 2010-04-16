// $Id$

#include <string>
#include "nuto/math/Matrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrixCSR<double>::GetTypeId() const
{
    return std::string("SparseMatrixCSRDouble");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrixCSR<int>::GetTypeId() const
{
    return std::string("SparseMatrixCSRInt");
}

template<>
void SparseMatrixCSR<int>::RemoveZeroEntries(double rAbsoluteTolerance, double rRelativeTolerance)
{

}
template<>
void SparseMatrixCSR<double>::RemoveZeroEntries(double rAbsoluteTolerance, double rRelativeTolerance)
{
    double tolerance =  rAbsoluteTolerance;
    if (rRelativeTolerance > 0)
    {
        double maxValue = 0;
        for (unsigned int entryCount = 0; entryCount < this->mValues.size(); entryCount++)
        {
            if (fabs(this->mValues[entryCount]) > maxValue)
            {
                maxValue = fabs(this->mValues[entryCount]);
            }
        }
        tolerance += rRelativeTolerance * maxValue;
    }

    int newPos = 0;
    int start = newPos;
    if (this->mOneBasedIndexing)
    {
        for (int row = 0; row < this->GetNumRows(); row++)
        {
            for (int pos = start; pos < this->mRowIndex[row + 1] - 1; pos++)
            {
                if (fabs(this->mValues[pos - 1]) > tolerance)
                {
                    this->mValues[newPos - 1] = this->mValues[pos];
                    this->mColumns[newPos - 1] = this->mColumns[pos];
                    newPos++;
                }
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
                if (fabs(this->mValues[pos]) > tolerance)
                {
                    this->mValues[newPos] = this->mValues[pos];
                    this->mColumns[newPos] = this->mColumns[pos];
                    newPos++;
                }
            }
            start = this->mRowIndex[row + 1];
            this->mRowIndex[row + 1] = newPos;
        }
    }
}

}
