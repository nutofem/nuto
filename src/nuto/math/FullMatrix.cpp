// $Id$

#include <boost/spirit/include/classic_core.hpp>
#include <fstream>
#include <iostream>

#include "nuto/math/dlapack.h"
#include "nuto/math/Matrix.h"
#include "nuto/math/FullMatrix.h"
#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>

namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string FullMatrix<double>::GetTypeId()const
{
    return std::string("FullMatrixDouble");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string FullMatrix<int>::GetTypeId()const
{
    return std::string("FullMatrixInt");
}

//! @brief ... converts any IntFullMatrix to an DoubleFullMatrix
//! @return    converted matrix
template<>
FullMatrix<double> FullMatrix<int>::Convert2double()
{
    FullMatrix<double> doubleMatrix(GetNumRows(),GetNumColumns());
    for (int count=0; count<GetNumColumns(); count++)
        for (int count2=0; count2<GetNumRows(); count2++)
            doubleMatrix(count2,count) = (double)mEigenMatrix(count2,count);

    return doubleMatrix;
}

//! @brief ... converts any DoubleFullMatrix to an DoubleFullMatrix
//! @return    converted matrix
template<>
FullMatrix<double> FullMatrix<double>::Convert2double()
{
    return Copy();
}

//! @brief ... converts any DoubleFullMatrix to an IntFullMatrix
//! @return    converted matrix
template<>
FullMatrix<int> FullMatrix<double>::Convert2int()
{
    FullMatrix<int> intMatrix(GetNumRows(),GetNumColumns());
    for (int count=0; count<GetNumColumns(); count++)
        for (int count2=0; count2<GetNumRows(); count2++)
            intMatrix(count2,count) = (int)mEigenMatrix(count2,count);

    return intMatrix;
}

//! @brief ... converts any IntFullMatrix to an IntFullMatrix
//! @return    converted matrix
template<>
FullMatrix<int> FullMatrix<int>::Convert2int()
{
    return Copy();
}

// import SLang matrix or vector
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
            this->mEigenMatrix(rowCount, colCount) = rowValues[colCount];
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

template<>
void FullMatrix<int>::ImportFromVtkASCIIFile(const char* rFileName)
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
        this->mEigenMatrix(count, 0) = imageValues[count];
    }
    // close file
   file.close();
}
template<>
void FullMatrix<double>::ImportFromVtkASCIIFile(const char* rfileName)
{
    throw MathException("[FullMatrix::importFromVtkASCIIFile] not implemented for double data-type.");
}

//! @brief elementwise absolute value of the matrix
template<>
FullMatrix<int> FullMatrix<int>::Abs() const
{
	return FullMatrix<int> ( mEigenMatrix.array().abs() );
}

template<>
 FullMatrix<double> FullMatrix<double>::Abs() const
{
	return FullMatrix<double> ( mEigenMatrix.array().abs() );
}
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
	FullMatrix<double> choleskyMatrix(this->mEigenMatrix);
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
	dpotrf_(&uplo, &dimMatrix, rInverse.mEigenMatrix.data(), &dimMatrix, &info);
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
	double* data = rInverse.mEigenMatrix.data();
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

//! @brief calculates the inverse of a matrix
//! @param rInverse ... inverse matrix
template<>
FullMatrix<int> FullMatrix<int>::Inverse() const
{
    throw MathException("[FullMatrix::Inverse] not implemented for integer data-type.");
}

//! @brief calculates the inverse of a matrix
//! @param rInverse ... inverse matrix
template<>
FullMatrix<double> FullMatrix<double>::Inverse() const
{
    return FullMatrix<double>(mEigenMatrix.inverse());
}


//! @brief calculates the inverse of a matrix
//! @param rInverse ... norm of matrix
template<>
double FullMatrix<int>::Norm() const
{
    throw MathException("[FullMatrix::Norm] not implemented for integer data-type.");
}

//! @brief ... calculates the norm of this matrix, i.e. for vectors the Euclidean norm
//! @return norm of this matrix
template<>
double FullMatrix<double>::Norm() const
{
    return mEigenMatrix.norm();
}


//! @brief calculates the eigenvalues
//! @param rEigenValues ... eigenvalues
template<>
void FullMatrix<int>::EigenValuesSymmetric(FullMatrix<double>& rEigenValues) const
{
    throw MathException("[FullMatrix::EigenValues] not implemented for integer data-type.");
}

//! @brief calculates the eigenvalues of a symmetric matrix
//! attention, of the matrix is not self adjoint (symmetric for real matrices), the result is wrong
//! @param rEigenValues ... eigenvalues
template<>
void FullMatrix<double>::EigenValuesSymmetric(FullMatrix<double>& rEigenValues) const
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> mySolver(mEigenMatrix,false);
    rEigenValues = FullMatrix<double>(mySolver.eigenvalues());
}

//! @brief calculates the eigenvalues
//! @param rEigenValues ... eigenvalues
template<>
void FullMatrix<int>::EigenVectorsSymmetric(FullMatrix<double>& rEigenVectors) const
{
    throw MathException("[FullMatrix::EigenValues] not implemented for integer data-type.");
}

//! @brief calculates the eigenvalues
//! @param rEigenValues ... the computed eigen vectors as a matrix of column vectors
template<>
void FullMatrix<double>::EigenVectorsSymmetric(FullMatrix<double>& rEigenVectors) const
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> mySolver(mEigenMatrix);
    rEigenVectors = FullMatrix<double>(mySolver.eigenvectors());
}

#ifdef ENABLE_SERIALIZATION
template<typename T>
void FullMatrix<T>::Save ( const std::string &filename, std::string rType)const
{
    try
    {
	//transform to uppercase
	std::transform(rType.begin(), rType.end(), rType.begin(), toupper);

	// open file
	std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
	if (!ofs.is_open())
	{
	    throw MathException ( std::string ( "[NuTo::FullMatrix::Save] error opening file." ) );
	}

	// write data
	std::string typeIdString ( GetTypeId() );
	if (rType=="BINARY")
	{
	    boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
	    oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="XML")
	{
	    boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
	    oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="TEXT")
	{
	    boost::archive::text_oarchive ota ( ofs, std::ios::binary );
	    ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else
	{
	    throw MathException ( "[NuTo::FullMatrix::Save]File type not implemented." );
	}

	// close file
	ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
	std::string s ( std::string ( "[NuTo::FullMatrix::Save]File save exception in boost - " ) +std::string ( e.what() ) );
	std::cout << s << "\n";
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
    catch ( ... )
    {
	throw MathException ( "[NuTo::FullMatrix::Save] Unhandled exception." );
    }
}

template<typename T>
void FullMatrix<T>::Restore ( const std::string &filename,  std::string rType)
{
    try
    {
	//transform to uppercase
	std::transform(rType.begin(), rType.end(), rType.begin(), toupper);

	// open file
	std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
	if (!ifs.is_open())
	{
	    throw MathException ( "[NuTo::FullMatrix::Restore] error opening file");
	}

	// read date
	std::string typeIdString;
	if (rType=="BINARY")
	{
	    boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
	    oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    if ( typeIdString!=GetTypeId() )
	    {
		throw MathException ( "[NuTo::FullMatrix::Restore]Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
	    }
	    oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="XML")
	{
	    boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
	    oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    if ( typeIdString!=GetTypeId() )
	    {
		throw MathException ( "[NuTo::FullMatrix::Restore]Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
	    }
	    oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else if (rType=="TEXT")
	{
	    boost::archive::text_iarchive ota ( ifs, std::ios::binary );
	    ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
	    if ( typeIdString!=GetTypeId() )
	    {
		throw MathException ( "[NuTo::FullMatrix::Restore]Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
	    }
	    ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
	}
	else
	{
	    throw MathException ( "[NuTo::FullMatrix::Restore]File type not implemented" );
	}

	// close file
	ifs.close();
    }
    catch ( MathException &e )
    {
	throw e;
    }
    catch ( std::exception &e )
    {
	throw MathException ( e.what() );
    }
    catch ( ... )
    {
	throw MathException ( "[NuTo::FullMatrix::Restore] Unhandled exception." );
    }
}


template void FullMatrix<int>::Save (const std::string&, std::string) const;
template void FullMatrix<int>::Restore (const std::string&, std::string);
template void FullMatrix<double>::Save (const std::string&, std::string) const;
template void FullMatrix<double>::Restore (const std::string&, std::string);

#endif // ENABLE_SERIALIZATION

}
