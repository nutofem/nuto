#include <boost/spirit/include/classic_core.hpp>
#include <fstream>
#include <iostream>
#include <string>

#include "math/SparseMatrix.h"
#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/MathException.h"

namespace NuTo
{

template<>
SparseMatrixCSRGeneral<int>::SparseMatrixCSRGeneral(const Eigen::MatrixXi& rFullMatrix, double rAbsoluteTolerance, double rRelativeTolerance): SparseMatrixCSR<int>(0,0)
{
    throw MathException("[SparseMatrixCSRGeneral::SparseMatrixCSRGeneral] conversion from full matrix not implemented for integers.");
    mNumColumns = 0;
}

template<>
SparseMatrixCSRGeneral<double>::SparseMatrixCSRGeneral(const Eigen::MatrixXd& rFullMatrix, double rAbsoluteTolerance, double rRelativeTolerance): SparseMatrixCSR<double>(0,0)
{
    this->Resize(rFullMatrix.rows(),rFullMatrix.cols());
    double tolerance = rAbsoluteTolerance;
    if (rRelativeTolerance > 1e-14)
    {
        double maxValue = 0;
        const double* values = rFullMatrix.data();
        for (int count = 0; count < rFullMatrix.rows() * rFullMatrix.cols(); count++)
        {
            if (std::abs(values[count]) > maxValue)
            {
                maxValue = std::abs(values[count]);
            }
        }
        tolerance += rRelativeTolerance * maxValue;
    }
    assert(this->mOneBasedIndexing == false);
    this->mRowIndex[0] = 0;
    for (int row = 0; row < rFullMatrix.rows(); row++)
    {
        for (int col = 0; col < rFullMatrix.cols(); col++)
        {
            if (std::abs(rFullMatrix(row,col)) > tolerance)
            {
                this->mValues.push_back(rFullMatrix(row,col));
                this->mColumns.push_back(col);
            }
        }
        assert(this->mValues.size() == this->mColumns.size());
        this->mRowIndex[row + 1] = this->mValues.size();
    }
    mNumColumns = rFullMatrix.cols();
}

template<>
void SparseMatrixCSRGeneral<int>::Gauss(Eigen::MatrixXi& rRhs, std::vector<int>& rMappingToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
{
    throw MathException("[SparseMatrixCSRGeneral::Gauss] not implemented for this data-type.");
}

template<>
void SparseMatrixCSRGeneral<int>::Gauss(SparseMatrixCSRGeneral<int>& rRhs, std::vector<int>& rMappingToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
{
    throw MathException("[SparseMatrixCSRGeneral::Gauss] not implemented for this data-type.");
}

template<>
void SparseMatrixCSRGeneral<double>::Gauss(Eigen::MatrixXd& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
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
    for (double value : mValues)
    {
        if (std::abs(value) > tolerance)
        {
            tolerance = std::abs(value);
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
            if (std::abs(this->mValues[pos]) > std::abs(pivot))
            {
                pivot = this->mValues[pos];
                swapCol = rMappingInitialToNewOrdering[this->mColumns[pos]];
            }
        }
        if (std::abs(pivot) < tolerance)
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
        for (int rhsCount = 0; rhsCount < rRhs.cols(); rhsCount++)
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
                    for (int rhsCount = 0; rhsCount < rRhs.cols(); rhsCount++)
                    {
                        rRhs(tmpRow, rhsCount) *= factor;
                    }

                    // add row to tmpRow
                    for (int rowPos = this->mRowIndex[row]; rowPos < this->mRowIndex[row + 1]; rowPos++)
                    {
                        this->AddValue(tmpRow, this->mColumns[rowPos], this->mValues[rowPos]);
                    }
                    for (int rhsCount = 0; rhsCount < rRhs.cols(); rhsCount++)
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
                if (std::abs(this->mValues[pos]) < tolerance)
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
        std::cout << "right-hand side after factorization \n" << rRhs << std::endl;
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
                    this->AddValue(row, this->mColumns[tmpPos], factor * this->mValues[tmpPos]);
                    //this is a dirty trick to check, if an element has been inserted in a previous row, which means, the position of the current value is
                    // shifted by one
                    if (oldsize != mValues.size())
                        tmpPos++;
                }
                for (int rhsCount = 0; rhsCount < rRhs.cols(); rhsCount++)
                {
                    rRhs(row, rhsCount) += factor * rRhs(column, rhsCount);
                }
            }
        }
        // remove zero entries
        for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
        {
        	if (std::abs(this->mValues[pos]) < tolerance)
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
        std::cout << "right-hand side after back substitution \n" << rRhs << std::endl;
        std::cout << "### end back substitution ###" << std::endl;
    }

    // renumber columns
    this->ReorderColumns(rMappingInitialToNewOrdering);
}


template<>
void SparseMatrixCSRGeneral<double>::Gauss(SparseMatrixCSRGeneral<double>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
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
    for (double value : mValues)
    {
        if (std::abs(value) > tolerance)
        {
            tolerance = std::abs(value);
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
            if (std::abs(this->mValues[pos]) > std::abs(pivot))
            {
                pivot = this->mValues[pos];
                swapCol = rMappingInitialToNewOrdering[this->mColumns[pos]];
            }
        }
        if (std::abs(pivot) < tolerance)
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

        for (int pos = rRhs.mRowIndex[row]; pos < rRhs.mRowIndex[row + 1]; pos++)
        {
        	rRhs.mValues[pos] *= invPivot;
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
                    for (int pos = rRhs.mRowIndex[tmpRow]; pos < rRhs.mRowIndex[tmpRow + 1]; pos++)
                    {
                    	rRhs.mValues[pos] *= factor;
                    }

                    // add row to tmpRow
                    for (int rowPos = this->mRowIndex[row]; rowPos < this->mRowIndex[row + 1]; rowPos++)
                    {
                        this->AddValue(tmpRow, this->mColumns[rowPos], this->mValues[rowPos]);
                    }
                    for (int rowPos = rRhs.mRowIndex[row]; rowPos < rRhs.mRowIndex[row + 1]; rowPos++)
                    {
                    	rRhs.AddValue(tmpRow, rRhs.mColumns[rowPos], rRhs.mValues[rowPos]);
                    }
                }
            }
        }

        // remove zero entries
        for (int tmpRow = row + 1; tmpRow < this->GetNumRows(); tmpRow++)
        {
            for (int pos = this->mRowIndex[tmpRow]; pos < this->mRowIndex[tmpRow + 1]; pos++)
            {
                if (std::abs(this->mValues[pos]) < tolerance)
                {
                    this->RemoveEntry(tmpRow, this->mColumns[pos]);
                    pos--;
                }
            }
        }
        // remove zero entries
        for (int tmpRow = row + 1; tmpRow < rRhs.GetNumRows(); tmpRow++)
        {
            for (int pos = rRhs.mRowIndex[tmpRow]; pos < rRhs.mRowIndex[tmpRow + 1]; pos++)
            {
                if (std::abs(rRhs.mValues[pos]) < tolerance)
                {
                	rRhs.RemoveEntry(tmpRow, rRhs.mColumns[pos]);
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
                    unsigned int oldsize = this->mValues.size();
                    this->AddValue(row, this->mColumns[tmpPos], factor * this->mValues[tmpPos]);
                    //this is a dirty trick to check, if an element has been inserted in a previous row, which means, the position of the current value is
                    // shifted by one
                    if (oldsize != this->mValues.size())
                        tmpPos++;
                }
                for (int tmpPos = rRhs.mRowIndex[column]; tmpPos < rRhs.mRowIndex[column + 1]; tmpPos++)
                {
                    unsigned int oldsize = rRhs.mValues.size();
                    rRhs.AddValue(row, rRhs.mColumns[tmpPos], factor * rRhs.mValues[tmpPos]);
                    //this is a dirty trick to check, if an element has been inserted in a previous row, which means, the position of the current value is
                    // shifted by one
                    if (oldsize != rRhs.mValues.size())
                        tmpPos++;
                }
            }
        }
        // remove zero entries
        for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
        {
        	if (std::abs(this->mValues[pos]) < tolerance)
            {
                this->RemoveEntry(row, this->mColumns[pos]);
                pos--;
            }
        }
        // remove zero entries
        for (int pos = rRhs.mRowIndex[row]; pos < rRhs.mRowIndex[row + 1]; pos++)
        {
        	if (std::abs(rRhs.mValues[pos]) < tolerance)
            {
        		rRhs.RemoveEntry(row, rRhs.mColumns[pos]);
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

template<>
void SparseMatrixCSRGeneral<int>::GetMaximumEigenvalueAndEigenvector(Eigen::VectorXi &rStart, int &maximumEigenvalue, double tol)
{
	throw MathException(__PRETTY_FUNCTION__, "not implemented for this data-type.");
}

template<>
void SparseMatrixCSRGeneral<double>::GetMaximumEigenvalueAndEigenvector(Eigen::VectorXd &rStart, double &maximumEigenvalue, double tol)
{
	int numRows = this->GetNumRows();

	assert(numRows != this->GetNumColumns());

	if(rStart.rows() != numRows)
	{
		rStart.resize(numRows);
		rStart.fill(1./std::sqrt(numRows + 1));
		rStart(numRows - 1) = std::sqrt(1 - (numRows - 1)/(numRows + 1));
	}

    Eigen::VectorXd y_k_1_star(rStart.rows());
    Eigen::VectorXd y_k_1(rStart.rows());
	double lambda_k_1 = 0., lambda_k_2 = 0.;
	double error = 0.;
	int i = 1;
	while(error > tol && i > 3)
	{
		y_k_1_star = this->operator*(rStart);
		maximumEigenvalue = rStart.dot(y_k_1_star);
		y_k_1 = y_k_1_star*(1./y_k_1_star.norm());

		rStart = y_k_1;

		if(i > 3)
		{
			double qk = (maximumEigenvalue - lambda_k_1)/(lambda_k_1 - lambda_k_2);
			double temp = std::abs((qk/(1-qk))*(maximumEigenvalue - lambda_k_1));
			error = temp/(temp+maximumEigenvalue);
		}

		lambda_k_2 = lambda_k_1;
		lambda_k_1 = maximumEigenvalue;
		i++;
	}
}

#ifdef ENABLE_SERIALIZATION
template <class T>
void SparseMatrixCSRGeneral<T>::Save ( const std::string &filename, std::string rType)const
{
    try
    {
        // transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int (*)(int))toupper);

        // open file
        std::ofstream ofs(filename.c_str(), std::ios_base::binary);
        if (!ofs.is_open())
        {
            throw MathException(__PRETTY_FUNCTION__, "Error opening file.");
        }

        // write data to file
        std::string typeIdString(this->GetTypeId());
        if (rType == "BINARY")
        {
            boost::archive::binary_oarchive oba(ofs, std::ios::binary);
            oba& boost::serialization::make_nvp("Object_type", typeIdString);
            oba& boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType == "XML")
        {
            boost::archive::xml_oarchive oxa(ofs, std::ios::binary);
            std::string tmpString(this->GetTypeId());
            oxa& boost::serialization::make_nvp("Object_type", typeIdString);
            oxa& boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType == "TEXT")
        {
            boost::archive::text_oarchive ota(ofs, std::ios::binary);
            ota& boost::serialization::make_nvp("Object_type", typeIdString);
            ota& boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MathException(__PRETTY_FUNCTION__, "File type not implemented.");
        }

        // close file
        ofs.close();
    }
    catch (boost::archive::archive_exception& e)
    {
        std::string s(__PRETTY_FUNCTION__ + "File save exception in boost - " + e.what());
        throw MathException(s);
    }
    catch (MathException& e)
    {
        throw;
    }
    catch (std::exception& e)
    {
        throw MathException(e.what());
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
	catch ( boost::archive::archive_exception& e )
	{
		std::string s ( std::string ( "[NuTo::SparseMatrixCSRGeneral::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
		throw MathException ( s );
	}
	catch ( MathException &e )
	{
        throw;
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
