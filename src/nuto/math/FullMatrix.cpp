// $Id$

#include <boost/spirit/include/classic_core.hpp>
#include <fstream>
#include <iostream>

#include "nuto/math/FullMatrix.h"
#include "nuto/math/dlapack.h"
#include "nuto/math/Matrix.h"

#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>

namespace NuTo
{
/*// import SLang matrix or vector
template<>
void FullMatrix<double>::ImportFromSLangText(const char* fileName)
{
    using namespace boost::spirit::classic;

    // open file
    std::ifstream file(fileName, std::ios::in);
    if (file.is_open() == false)
    {
        throw MathException("[FullMatrix::importFromSLang]error opening file.");
    }

    // read header
    unsigned int SLangVersion(0), objectType(SLANG_NOTYPE), objectKind(SLANG_NOKIND), objectNumRows(0), objectNumColumns(0), objectNumEntries(0);
    this->ImportFromSLangTextReadHeader(file, SLangVersion, objectType, objectKind, objectNumRows, objectNumColumns, objectNumEntries);

    if (SLangVersion != 511)
    {
        throw MathException("[FullMatrix::importFromSLang]unsupported SLang version.");
    }

    if (objectType != SLANG_REAL)
    {
        throw MathException("[FullMatrix::importFromSLang]object data must be of type DOUBLE.");
    }
    if (objectKind != SLANG_VECTOR && objectKind != SLANG_MATRIX)
    {
        throw MathException("[FullMatrix::importFromSLang]object must be a MATRIX or a VECTOR.");
    }
    this->Resize(objectNumRows, objectNumColumns);

    // read entries
    std::vector<double> rowValues;
    rowValues.reserve(objectNumColumns);
    for (unsigned int rowCount = 0; rowCount < objectNumRows; rowCount++)
    {
        std::string line;
        getline (file, line);
        rowValues.clear();
        if (parse(line.c_str(),(*space_p >> *(real_p[push_back_a(rowValues)] >> *space_p))).full == false)
        {
            throw MathException("[FullMatrix::importFromSLang]error reading entries.");
        }
        if (rowValues.size() != objectNumColumns)
        {
            throw MathException("[FullMatrix::importFromSLang]invalid number of values in row.");
        }
        for (unsigned int colCount = 0; colCount < objectNumColumns; colCount++)
        {
            (*this)(rowCount, colCount) = rowValues[colCount];
        }
    }

    // close file
    file.close();
}

template<>
void FullMatrix<int>::ImportFromSLangText(const char* fileName)
{
    throw MathException("[FullMatrix::importFromSLang] not implemented for this data-type.");
}
*/

/*
template<>
void FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic>::ImportFromVtkASCIIFile(const char* rFileName)
{
    using namespace boost::spirit::classic;

    // open file
    std::ifstream file(rFileName, std::ios::in);
    if (file.is_open() == false)
    {
        throw MathException("[FullMatrix::ImportFromVtkASCIIFile] error opening file.");
    }
    //skip read header, is done in StructureGrid.cpp
    // read first four lines
    int numEntries(0);

    std::string line;
    for (int count=0;count<7;count++)
    {
        getline (file, line);
    }

    // read number of entries
    getline (file, line);
    if (parse(line.c_str(),("POINT_DATA ">> int_p[assign_a(numEntries)] >>  *space_p)).full == false)
           {
               throw MathException("[Matrix::importFromVtkASCIIFile]error reading number of entries.");
           }
    // read data type
    getline (file, line);
    // read empty line
    getline (file, line);

  // read entries
    std::vector<int> imageValues;
    int value;
    imageValues.reserve(numEntries);

    while(getline(file,line))
    {
    	std::istringstream iss(line);
    	while(iss >> value)
    	{
			imageValues.push_back(value);
    	}
    }
    for (int count = 0; count<numEntries; count++)
    {
        (*this)(count, 0) = imageValues[count];
    }
    // close file
   file.close();
}
template<>
void FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>::ImportFromVtkASCIIFile(const char* rfileName)
{
    throw MathException("[FullMatrix::importFromVtkASCIIFile] not implemented for double data-type.");
}

*/

/*
template<>
void FullMatrix<int>::SolveCholeskyLapack(const FullMatrix<double>& rRHS, FullMatrix<double>& rSolution) const
{
	throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] not implemented for integer data-type.");
}

template<>
void FullMatrix<double>::SolveCholeskyLapack(const FullMatrix<double>& rRHS, FullMatrix<double>& rSolution) const
{
#ifdef ENABLE_MKL
	// check matrix
	if(this->GetNumColumns() != this->GetNumRows())
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid shape of the coefficient matrix.");
	}
	if(this->GetNumColumns() < 1)
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid dimension of the coefficient matrix.");
	}
	// check right-hand-side vectors
	if(this->GetNumRows() != rRHS.GetNumRows())
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid number of rows in right-hand-side matrix.");
	}
	if(rRHS.GetNumColumns() < 1)
	{
		throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] invalid number of columns in right-hand-side matrix.");
	}

	// copy coefficient matrix and perform cholesky factorization
	FullMatrix<double> choleskyMatrix(*this);
	int dimCholeskyMatrix = choleskyMatrix.GetNumColumns();
	char uplo('L');
	int info(0);
	dpotrf_(&uplo, &dimCholeskyMatrix, choleskyMatrix.mEigenMatrix.data(), &dimCholeskyMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::SolveCholeskyLapack] Error calculating Cholesky factorization.");
		if(info < 1)
		{
			info *= -1;
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		}
		else
		{
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The leading minor of order "+infoStream.str()+"is not positive definite.";
		}
		throw MathException(errorMessage);
	}

	// copy right-hand-side matrix and solve linear system of equations using factorized choolesky matrix
	rSolution = rRHS;
	int numRHSVectors = rSolution.GetNumColumns();
	dpotrs_(&uplo, &dimCholeskyMatrix, &numRHSVectors, choleskyMatrix.mEigenMatrix.data(), &dimCholeskyMatrix, rSolution.mEigenMatrix.data(), &dimCholeskyMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::SolveCholeskyLapack] Error solving linear system of equations.");
		std::stringstream infoStream;
		infoStream << -1*info;
		errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		throw MathException(errorMessage);
	}
#else //ENABLE_MKL
	throw MathException("[NuTo::FullMatrix::SolveCholeskyLapack] lapack package not enabled.");
#endif //ENABLE_MKL
}

template<>
void FullMatrix<int>::InverseCholeskyLapack(FullMatrix<double>& rInverse) const
{
	throw MathException("[FullMatrix::InverseCholeskyLapack] not implemented for integer data-type.");
}

template<>
void FullMatrix<double>::InverseCholeskyLapack(FullMatrix<double>& rInverse) const
{
#ifdef ENABLE_MKL
	// check matrix
	if(this->GetNumColumns() != this->GetNumRows())
	{
		throw MathException("[NuTo::FullMatrix::InverseCholeskyLapack] invalid shape of the coefficient matrix.");
	}
	if(this->GetNumColumns() < 1)
	{
		throw MathException("[NuTo::FullMatrix::InverseCholeskyLapack] invalid dimension of the coefficient matrix.");
	}

	// copy matrix and perform cholesky factorization
	rInverse = *this;
	int dimMatrix = rInverse.GetNumColumns();
	char uplo('L');
	int info(0);
	dpotrf_(&uplo, &dimMatrix, rInverse.data(), &dimMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::InverseCholeskyLapack] Error calculating Cholesky factorization.");
		if(info < 1)
		{
			info *= -1;
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		}
		else
		{
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The leading minor of order "+infoStream.str()+"is not positive definite.";
		}
		throw MathException(errorMessage);
	}

	// calculate inverse matrix
	dpotri_(&uplo, &dimMatrix, rInverse.mEigenMatrix.data(), &dimMatrix, &info);
	if(info != 0)
	{
		std::string errorMessage("[NuTo::FullMatrix::InverseCholeskyLapack] Error calculating inverse matrix.");
		if(info < 1)
		{
			info *= -1;
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The "+infoStream.str()+"-th argument had an illegal value.";
		}
		else
		{
			std::stringstream infoStream;
			infoStream << info;
			errorMessage += " The ("+infoStream.str()+","+infoStream.str()+")-th element of the factor U or L is zero, and the inverse could not be computed.";
		}
		throw MathException(errorMessage);
	}

	// copy values from lower triangle to upper triangle
	double* data = rInverse.data();
	for(int col = 1; col < dimMatrix; col++)
	{
		for(int row = 0; row < col; row++)
		{
			data[col * dimMatrix + row] = data[row * dimMatrix + col];
		}
	}
#else //ENABLE_MKL
	throw MathException("[NuTo::FullMatrix::InverseCholeskyLapack] lapack package not enabled.");
#endif //ENABLE_MKL
}
*/

template<>
std::string GetTypeIdBaseDataType<int>()
{
	return std::string("Int");
}

template<>
std::string GetTypeIdBaseDataType<double>()
{
	return std::string("Double");
}



#ifdef ENABLE_SERIALIZATION

template void FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>::Save (const std::string&, std::string) const;
template void FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>::Restore (const std::string&, std::string);
template void FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>::Save (const std::string&, std::string) const;
template void FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>::Restore (const std::string&, std::string);

#endif // ENABLE_SERIALIZATION
} //NAMESPACE NUTO


