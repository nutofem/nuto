// $Id: SparseMatrixCSRVector2General.h 235 2010-04-22 09:25:38Z arnold2 $

#ifndef SPARSE_MATRIX_CSR_VECTOR2_GENERAL_H
#define SPARSE_MATRIX_CSR_VECTOR2_GENERAL_H

#include "nuto/math/Matrix.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... class for general sparse matrices which are stored in CSR format
template <class T>
class SparseMatrixCSRVector2General : public SparseMatrixCSRVector2<T>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns_ ... number of columns
    SparseMatrixCSRVector2General(int rNumRows_=0, int rNumColumns_=0) : SparseMatrixCSRVector2<T>(rNumRows_,rNumColumns_)
    {
    }

    //! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
    //! @param rFullMatrix ... input matrix (full storage)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
    SparseMatrixCSRVector2General(FullMatrix<T>& rFullMatrix, double rAbsoluteTolerance = 0, double rRelativeTolerance = 1e-14): SparseMatrixCSRVector2<T>(rFullMatrix.GetNumRows(),rFullMatrix.GetNumColumns())
    {
        double tolerance = rAbsoluteTolerance;
        if (rRelativeTolerance > 1e-14)
        {
            T maxValue(0), minValue(0);
            this->Max(maxValue);
            this->Min(minValue);
            if (fabs(maxValue)>fabs(minValue))
                tolerance += rRelativeTolerance * fabs(maxValue);
            else
                tolerance += rRelativeTolerance * fabs(minValue);

        }

        for (int row = 0; row < rFullMatrix.GetNumRows(); row++)
        {
            for (int col = 0; col < rFullMatrix.GetNumColumns(); col++)
            {
                if (rFullMatrix(row,col) > tolerance)
                {
                    this->AddEntry(rFullMatrix(row,col),row,col);
                }
            }
        }
    }

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const;

    //! @brief ... returns whether the matrix is symmetric or unsymmetric
    //! @return true if the matrix is symmetric and false if the matrix is unsymmetric
    bool IsSymmetric() const
    {
        return false;
    }

    //! @brief ... add nonzero entry to matrix
    //! @param rRow ... row of the nonzero entry (zero based indexing!!!)
    //! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
    //! @param rValue ... value of the nonzero entry
    void AddEntry(int rRow, int rColumn, T rValue)
    {
        // check for overflow
        assert(rRow < INT_MAX);
        assert(rColumn < INT_MAX);
        assert(rRow >= 0);
        assert(rColumn >= 0);

        // check bounds
        if (rRow > (int)this->mValues.size() || rRow<0)
        {
            throw MathException("[SparseMatrixCSRVector2General::addEntry] row index is out of bounds.");
        }
        if (rColumn >= this->mNumColumns || rColumn<0)
        {
            throw MathException("[SparseMatrixCSRVector2General::addEntry] column index is out of bounds.");
        }
        if (this->mOneBasedIndexing)
        {
        	rColumn++;
        }
		// find position in matrix
		bool existingValue(false);
		int col_pos;
		for (col_pos = 0; (col_pos < (int)this->mColumns[rRow].size()); col_pos++)
		{
			if (this->mColumns[rRow][col_pos] < rColumn)
				continue;
			else
			{
				if (this->mColumns[rRow][col_pos] == rColumn)
				{
					existingValue = true;
				}
				break;
			}
		}

		// add value
		if (existingValue)
		{
			// add to existing value
			this->mValues[rRow][col_pos] += rValue;
		}
		else
		{
			// insert new value
			this->mColumns[rRow].insert(this->mColumns[rRow].begin() + col_pos, rColumn);
			this->mValues[rRow].insert(this->mValues[rRow].begin() + col_pos, rValue);
		}
    }
    //! @brief ... import matrix from slang object stored in  a text file
    //! @param rFileName ... file name
    void ImportFromSLangText(const char* rFileName)
    {
    	throw MathException("NuTo::SparseMatrixCSRVector2General::ImportFromSLangText] to be implemented.");
    }

    //! @brief ... write nonzero matrix entries into a full matrix
    //! @param rFullMatrix ... the full matrix
    void WriteEntriesToFullMatrix(FullMatrix<T>& rFullMatrix) const
    {
        if (this->mOneBasedIndexing)
        {
            for (unsigned int row=0; row<this->mColumns.size(); row++)
            {
            	for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
            	{
            		rFullMatrix(row, this->mColumns[row][col_count]-1) = this->mValues[row][col_count];
            	}
            }
        }
        else
        {
            for (unsigned int row=0; row<this->mValues.size(); row++)
            {
            	for (unsigned int col_count=0; col_count<this->mColumns.size(); col_count++)
            	{
            		rFullMatrix(row, this->mColumns[row][col_count]) = this->mValues[row][col_count];
            	}
            }
        }
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const
    {
        try
         {
             //transform to uppercase
             std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

             // open file
             std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
             if(! ofs.is_open())
             {
                 throw MathException("[NuTo::SparseMatrixCSRVector2General::Save] Error opening file.");
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
                 throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Save] File type not implemented." );
             }

             // close file
             ofs.close();
         }
         catch ( boost::archive::archive_exception e )
         {
             std::string s ( std::string ( "[NuTo::SparseMatrixCSRVector2General::Save]File save exception in boost - " ) + std::string ( e.what() ) );
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

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType)
    {
        try
        {
            //transform to uppercase
            std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

            // open file
            std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
            if(! ifs.is_open())
            {
                throw MathException("[NuTo::SparseMatrixCSRVector2General::Restore] Error opening file.");
            }

            std::string typeIdString;
            if (rType=="BINARY")
            {
                boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
                oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
                if ( typeIdString != this->GetTypeId() )
                {
                    throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
                }
                oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
            }
            else if (rType=="XML")
            {
                boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
                oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
                if ( typeIdString != this->GetTypeId() )
                {
                    throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
                }
                oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
            }
            else if (rType=="TEXT")
            {
                boost::archive::text_iarchive ota ( ifs, std::ios::binary );
                ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
                if ( typeIdString != this->GetTypeId() )
                {
                    throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
                }
                ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
            }
            else
            {
                throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore]File type not implemented" );
            }
            // close file
            ifs.close();
        }
        catch ( boost::archive::archive_exception e )
        {
            std::string s ( std::string ( "[NuTo::SparseMatrixCSRVector2General::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
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

#endif // ENABLE_SERIALIZATION
protected:
};
}
#endif // SPARSE_MATRIX_CSR_VECTOR2_GENERAL_H
