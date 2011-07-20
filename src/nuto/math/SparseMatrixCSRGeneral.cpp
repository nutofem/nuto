// $Id$

#include <boost/spirit/include/classic_core.hpp>
#include <fstream>
#include <iostream>
#include <string>

#include "nuto/math/Matrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/MathException.h"

namespace NuTo
{

template<>
SparseMatrixCSRGeneral<int>::SparseMatrixCSRGeneral(const FullMatrix<int>& rFullMatrix, double rAbsoluteTolerance, double rRelativeTolerance): SparseMatrixCSR<int>(0,0)
{
    throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] conversion from full matrix not implemented for integers.");
}

template<>
SparseMatrixCSRGeneral<double>::SparseMatrixCSRGeneral(const FullMatrix<double>& rFullMatrix, double rAbsoluteTolerance, double rRelativeTolerance): SparseMatrixCSR<double>(0,0)
{
    this->Resize(rFullMatrix.GetNumRows(),rFullMatrix.GetNumColumns());
    double tolerance = rAbsoluteTolerance;
    if (rRelativeTolerance > 1e-14)
    {
        double maxValue = 0;
        const double* values = rFullMatrix.GetEigenMatrix().data();
        for (int count = 0; count < rFullMatrix.GetNumRows() * rFullMatrix.GetNumColumns(); count++)
        {
            if (fabs(values[count]) > maxValue)
            {
                maxValue = fabs(values[count]);
            }
        }
        tolerance += rRelativeTolerance * maxValue;
    }
    assert(this->mOneBasedIndexing == false);
    this->mRowIndex[0] = 0;
    for (int row = 0; row < rFullMatrix.GetNumRows(); row++)
    {
        for (int col = 0; col < rFullMatrix.GetNumColumns(); col++)
        {
            if (rFullMatrix(row,col) > tolerance)
            {
                this->mValues.push_back(rFullMatrix(row,col));
                this->mColumns.push_back(col);

            }
        }
        assert(this->mValues.size() == this->mColumns.size());
        this->mRowIndex[row + 1] = this->mValues.size();
    }
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrixCSRGeneral<double>::GetTypeId()const
{
    return std::string("SparseMatrixCSRGeneralDouble");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrixCSRGeneral<int>::GetTypeId()const
{
    return std::string("SparseMatrixCSRGeneralInt");
}

template<>
void SparseMatrixCSRGeneral<double>::ImportFromSLangText(const char* rFileName)
{
    using namespace boost::spirit::classic;

    // open file
    std::ifstream file(rFileName, std::ios::in);
    if (file.is_open() == false)
    {
        throw MathException("[SparseMatrixCSRGeneral::importFromSLangText]error opening file.");
    }

    // read header
    unsigned int SLangVersion(0), objectType(SLANG_NOTYPE), objectKind(SLANG_NOKIND), objectNumRows(0), objectNumColumns(0), objectNumEntries(0);
    this->ImportFromSLangTextReadHeader(file, SLangVersion, objectType, objectKind, objectNumRows, objectNumColumns, objectNumEntries);

    // check SLang version
    if (SLangVersion != 511)
    {
        throw MathException("[SparseMatrixCSRGeneral::importFromSLangText]unsupported SLang version.");
    }

    // check object type and kind
    if (objectType != SLANG_REAL)
    {
        throw MathException("[SparseMatrixCSRGeneral::importFromSLangText]object data must be of type DOUBLE.");
    }
    if (objectKind != SLANG_GENERAL_COMPACT_MATRIX)
    {
        throw MathException("[SparseMatrixCSRGeneral::importFromSLangText] SLang object must be a GENERAL_COMPACT_MATRIX.");
    }
    this->Resize(objectNumRows,objectNumColumns);
    this->Reserve(objectNumEntries);

    // switch to one based indexing
    bool oldOneBasedIndexing = this->mOneBasedIndexing;
    this->SetOneBasedIndexing();

    // read nonzero entries
    unsigned int row = 1;
    this->mRowIndex[0] = 1;
    for (unsigned int entryCount = 0; entryCount < objectNumEntries; entryCount++)
    {
        std::string line;
        getline (file, line);
        unsigned int curRow(0), curColumn(0);
        double curValue(0);
        if (parse(line.c_str(),(uint_p[assign_a(curRow)] >> real_p[assign_a(curColumn)] >> real_p[assign_a(curValue)]),space_p).full == false)
        {
            throw MathException("[SparseMatrixCSRGeneral::importFromSLangText]error reading nonzero matrix entries.");
        }
        if (curRow < row)
        {
            throw MathException("[SparseMatrixCSRGeneral::importFromSLangText]invalid sorting of compressed matrix.");
        }
        else
        {
            this->mColumns.push_back(curColumn);
            this->mValues.push_back(curValue);
            if (curRow > row)
            {
                while (row < curRow)
                {
                    this->mRowIndex[row] = entryCount + 1;
                    row++;
                }
            }
        }
    }
    while (row < this->mRowIndex.size())
    {
        this->mRowIndex[row] = objectNumEntries + 1;
        row++;
    }

    // set to old indexing
    if (!oldOneBasedIndexing)
    {
        this->SetZeroBasedIndexing();
    }

    // close file
    file.close();
}

template<>
void SparseMatrixCSRGeneral<int>::ImportFromSLangText(const char* rFileName)
{
    throw MathException("[SparseMatrixCSRGeneral::importFromSLang] not implemented for this data-type.");
}

template<>
void SparseMatrixCSRGeneral<int>::Gauss(FullMatrix<int>& rRhs, std::vector<int>& rMappingToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
{
    throw MathException("[SparseMatrixCSRGeneral::Gauss] not implemented for this data-type.");
}

template<>
void SparseMatrixCSRGeneral<double>::Gauss(FullMatrix<double>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
{
    // initialize help vectors for reordering
    rMappingNewToInitialOrdering.resize(this->GetNumColumns());
    rMappingInitialToNewOrdering.resize(this->GetNumColumns());
    for (int colCount = 0; colCount < this->GetNumColumns(); colCount++)
    {
        rMappingNewToInitialOrdering[colCount] = colCount;
        rMappingInitialToNewOrdering[colCount] = colCount;
    }

    // calculate tolerance
    double tolerance = 0;
    for (unsigned int valCount = 0; valCount < this->mValues.size(); valCount++)
    {
        if (fabs(this->mValues[valCount]) > tolerance)
        {
            tolerance = fabs(this->mValues[valCount]);
        }
    }
    tolerance *= rRelativeTolerance;

    if (this->mVerboseLevel > 0)
    {
        std::cout << "### start factorization ###" << std::endl;
    }
    for (int row = 0; row < this->GetNumRows(); row++)
    {
        // column pivoting
        // search for largest element in column
        double pivot = 0;
        int swapCol = 0;
        for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
        {
            if (fabs(this->mValues[pos]) > fabs(pivot))
            {
                pivot = this->mValues[pos];
                swapCol = rMappingInitialToNewOrdering[this->mColumns[pos]];
            }
        }
        if (fabs(pivot) < tolerance)
        {
            throw MathException("[SparseMatrixCSRGeneral<double>::Gauss] equation system is linear dependent.");
        }

        // now swap the columns
        int tmp = rMappingNewToInitialOrdering[row];
        rMappingNewToInitialOrdering[row] = rMappingNewToInitialOrdering[swapCol];
        rMappingNewToInitialOrdering[swapCol] = tmp;
        rMappingInitialToNewOrdering[rMappingNewToInitialOrdering[row]] = row;
        rMappingInitialToNewOrdering[rMappingNewToInitialOrdering[swapCol]] = swapCol;

        // unit value on the diagonal
        double invPivot = 1./pivot;
        for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
        {
            this->mValues[pos] *= invPivot;
        }
        for (int rhsCount = 0; rhsCount < rRhs.GetNumColumns(); rhsCount++)
        {
            rRhs(row, rhsCount) *= invPivot;
        }

        // linear combination of rows
        for (int tmpRow = row + 1; tmpRow < this->GetNumRows(); tmpRow++)
        {
            if (this->mRowIndex[tmpRow] == this->mRowIndex[tmpRow + 1])
            {
                throw MathException("[SparseMatrixCSRGeneral<double>::Gauss] equation system is linear dependent.");
            }
            else
            {
                int tmpPos = this->mRowIndex[tmpRow];
                for (; tmpPos < this->mRowIndex[tmpRow + 1]; tmpPos++)
                {
                    if (rMappingInitialToNewOrdering[this->mColumns[tmpPos]] == row)
                    {
                        break;
                    }
                }
                if (tmpPos < this->mRowIndex[tmpRow + 1])
                {
                    // scale tmpRow
                    double factor = -1./this->mValues[tmpPos];
                    for (int pos = this->mRowIndex[tmpRow]; pos < this->mRowIndex[tmpRow + 1]; pos++)
                    {
                        this->mValues[pos] *= factor;
                    }
                    for (int rhsCount = 0; rhsCount < rRhs.GetNumColumns(); rhsCount++)
                    {
                        rRhs(tmpRow, rhsCount) *= factor;
                    }

                    // add row to tmpRow
                    for (int rowPos = this->mRowIndex[row]; rowPos < this->mRowIndex[row + 1]; rowPos++)
                    {
                        this->AddEntry(tmpRow, this->mColumns[rowPos], this->mValues[rowPos]);
                    }
                    for (int rhsCount = 0; rhsCount < rRhs.GetNumColumns(); rhsCount++)
                    {
                        rRhs(tmpRow, rhsCount) += rRhs(row, rhsCount);
                    }
                }
            }
        }

        // remove zero entries
        for (int tmpRow = row + 1; tmpRow < this->GetNumRows(); tmpRow++)
        {
            for (int pos = this->mRowIndex[tmpRow]; pos < this->mRowIndex[tmpRow + 1]; pos++)
            {
                if (fabs(this->mValues[pos]) < tolerance)
                {
                    this->RemoveEntry(tmpRow, this->mColumns[pos]);
                    pos--;
                }
            }
        }
    }

    if (this->mVerboseLevel > 0)
    {
        std::cout << "column order after factorization" << std::endl;
        for (int colCount = 0; colCount < this->GetNumColumns(); colCount++)
        {
            std::cout << " " << rMappingNewToInitialOrdering[colCount];
        }
        std::cout << std::endl;
        std::cout << "matrix after factorization" << std::endl;
        for (int row = 0; row < this->GetNumRows(); row++)
        {
            for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row+1]; pos++)
            {
                std::cout << " row: " << row << " column: " << rMappingInitialToNewOrdering[this->mColumns[pos]] << " value: " << this->mValues[pos] << std::endl;
            }
        }
        std::cout << "right-hand side after factorization" << std::endl;
        rRhs.Info();
        std::cout << "### end factorization ###" << std::endl;
    }

    // back substitution (start with the last but one row
    if (this->mVerboseLevel > 0)
    {
        std::cout << "### start back substitution ###" << std::endl;
    }
    for (int row = this->GetNumRows() - 2; row >= 0; row--)
    {
        for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
        {
            int column = rMappingInitialToNewOrdering[this->mColumns[pos]];
            assert(column >= row);
            if ((row != column) && (column < this->GetNumRows()))
            {
                //linear combination of rows
                double factor = -this->mValues[pos];
                for (int tmpPos = this->mRowIndex[column]; tmpPos < this->mRowIndex[column + 1]; tmpPos++)
                {
                    unsigned int oldsize = mValues.size();
                    this->AddEntry(row, this->mColumns[tmpPos], factor * this->mValues[tmpPos]);
                    //this is a dirty trick to check, if an element has been inserted in a previous row, which means, the position of the current value is
                    // shifted by one
                    if (oldsize != mValues.size())
                        tmpPos++;
                }
                for (int rhsCount = 0; rhsCount < rRhs.GetNumColumns(); rhsCount++)
                {
                    rRhs(row, rhsCount) += factor * rRhs(column, rhsCount);
                }
            }
        }
        // remove zero entries
        for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
        {
        	if (fabs(this->mValues[pos]) < tolerance)
            {
                this->RemoveEntry(row, this->mColumns[pos]);
                pos--;
            }
        }
    }

    if (this->mVerboseLevel > 0)
    {
        std::cout << "matrix after back substitution" << std::endl;
        for (int row = 0; row < this->GetNumRows(); row++)
        {
            for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row+1]; pos++)
            {
                std::cout << " row: " << row << " column: " << rMappingInitialToNewOrdering[this->mColumns[pos]] << " value: " << this->mValues[pos] << std::endl;
            }
        }
        std::cout << "right-hand side after back substitution" << std::endl;
        rRhs.Info();
        std::cout << "### end back substitution ###" << std::endl;
    }

    // renumber columns
    this->ReorderColumns(rMappingInitialToNewOrdering);
}

#ifdef ENABLE_SERIALIZATION
template <class T>
void SparseMatrixCSRGeneral<T>::Save ( const std::string &filename, std::string rType)const
{
	try
	 {
		 //transform to uppercase
		 std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

		 // open file
		 std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
		 if(! ofs.is_open())
		 {
			 throw MathException("[NuTo::SparseMatrixCSRGeneral::Save] Error opening file.");
		 }

		 // write data to file
		 std::string typeIdString(this->GetTypeId());
		 if (rType=="BINARY")
		 {
			 boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
			 oba & boost::serialization::make_nvp ("Object_type", typeIdString );
			 oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		 }
		 else if (rType=="XML")
		 {
			 boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
			 std::string tmpString(this->GetTypeId());
			 oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
			 oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		 }
		 else if (rType=="TEXT")
		 {
			 boost::archive::text_oarchive ota ( ofs, std::ios::binary );
			 ota & boost::serialization::make_nvp("Object_type", typeIdString );
			 ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		 }
		 else
		 {
			 throw MathException ( "[NuTo::SparseMatrixCSRGeneral::Save] File type not implemented." );
		 }

		 // close file
		 ofs.close();
	 }
	 catch ( boost::archive::archive_exception e )
	 {
		 std::string s ( std::string ( "[NuTo::SparseMatrixCSRGeneral::Save]File save exception in boost - " ) + std::string ( e.what() ) );
		 throw MathException ( s );
	 }
	 catch ( MathException &e )
	 {
		 throw e;
	 }
	 catch ( std::exception &e )
	 {
		 throw MathException ( e.what() );
	 }
}

template <class T>
void SparseMatrixCSRGeneral<T>::Restore ( const std::string &filename,  std::string rType)
{
	try
	{
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

		// open file
		std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
		if(! ifs.is_open())
		{
			throw MathException("[NuTo::SparseMatrixCSRGeneral::Restore] Error opening file.");
		}

		std::string typeIdString;
		if (rType=="BINARY")
		{
			boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
			oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
			if ( typeIdString != this->GetTypeId() )
			{
				throw MathException ( "[NuTo::SparseMatrixCSRGeneral::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
			}
			oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		}
		else if (rType=="XML")
		{
			boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
			oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
			if ( typeIdString != this->GetTypeId() )
			{
				throw MathException ( "[NuTo::SparseMatrixCSRGeneral::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
			}
			oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		}
		else if (rType=="TEXT")
		{
			boost::archive::text_iarchive ota ( ifs, std::ios::binary );
			ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
			if ( typeIdString != this->GetTypeId() )
			{
				throw MathException ( "[NuTo::SparseMatrixCSRGeneral::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
			}
			ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		}
		else
		{
			throw MathException ( "[NuTo::SparseMatrixCSRGeneral::Restore]File type not implemented" );
		}
		// close file
		ifs.close();
	}
	catch ( boost::archive::archive_exception e )
	{
		std::string s ( std::string ( "[NuTo::SparseMatrixCSRGeneral::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
		throw MathException ( s );
	}
	catch ( MathException &e )
	{
		throw e;
	}
	catch ( std::exception &e )
	{
		throw MathException ( e.what() );
	}
}

template void SparseMatrixCSRGeneral<int>::Save (const std::string &, std::string) const;
template void SparseMatrixCSRGeneral<int>::Restore (const std::string &, std::string);
template void SparseMatrixCSRGeneral<double>::Save (const std::string &, std::string) const;
template void SparseMatrixCSRGeneral<double>::Restore (const std::string &, std::string);

#endif // ENABLE_SERIALIZATION

}
