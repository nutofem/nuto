// $Id: SparseMatrixCSRVector2.cpp 207 2009-12-18 08:08:29Z eckardt4 $

#include <string>
#include "nuto/math/Matrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrixCSRVector2<double>::GetTypeId() const
{
    return std::string("SparseMatrixCSRVector2Double");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrixCSRVector2<int>::GetTypeId() const
{
    return std::string("SparseMatrixCSRVector2Int");
}
template<>
int SparseMatrixCSRVector2<int>::RemoveZeroEntries(double rAbsoluteTolerance, double rRelativeTolerance)
{
	   double tolerance =  rAbsoluteTolerance;
	    if (rRelativeTolerance > 0)
	    {
	        int maxValue,minValue;
	        Max(maxValue);
	        Min(minValue);

	        if (fabs(maxValue)>fabs(minValue))
	            tolerance += rRelativeTolerance * maxValue;
	    }

	    int numRemoved(0);
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			int numRemovedPerRow;
			unsigned int newPos(0);
			std::vector<int>& thisValueVec(this->mValues[row]);
			std::vector<int>& thisColumnVec(this->mColumns[row]);
			for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
			{
				if (fabs(thisValueVec[pos]) > tolerance)
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

template<>
int SparseMatrixCSRVector2<double>::RemoveZeroEntries(double rAbsoluteTolerance, double rRelativeTolerance)
{
   double tolerance =  rAbsoluteTolerance;
	if (rRelativeTolerance > 0)
	{
		double maxValue,minValue;
		Max(maxValue);
		Min(minValue);

		if (fabs(maxValue)>fabs(minValue))
			tolerance += rRelativeTolerance * maxValue;
	}

	int numRemoved(0);
	for (int row = 0; row < this->GetNumRows(); row++)
	{
		unsigned int numRemovedPerRow(0);
		unsigned int newPos(0);
		std::vector<double>& thisValueVec(this->mValues[row]);
		std::vector<int>& thisColumnVec(this->mColumns[row]);
		for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
		{
			if (fabs(thisValueVec[pos]) > tolerance)
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

}
