// $Id$

#include <boost/spirit/include/classic_core.hpp>
#include <fstream>
#include <iostream>

#include "nuto/math/Matrix.h"
#include "nuto/math/FullMatrix.h"
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

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixShort
template<>
std::string FullMatrix<short>::GetTypeId()const
{
    return std::string("FullMatrixShort");
}

//! @brief ... converts any ShortFullMatrix to an DoubleFullMatrix
//! @return    converted matrix
template<>
FullMatrix<double> FullMatrix<short>::Convert2double()
{
    FullMatrix<double> doubleMatrix(GetNumRows(),GetNumColumns());
    for (int count=0; count<GetNumColumns(); count++)
        for (int count2=0; count2<GetNumRows(); count2++)
            doubleMatrix(count2,count) = (double)mEigenMatrix(count2,count);

    return doubleMatrix;
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

//! @brief ... converts any ShortFullMatrix to an IntFullMatrix
//! @return    converted matrix
template<>
FullMatrix<int> FullMatrix<short>::Convert2int()
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

//! @brief ... converts any DoubleFullMatrix to an ShortFullMatrix
//! @return    converted matrix
template<>
FullMatrix<short> FullMatrix<double>::Convert2short()
{
    FullMatrix<short> shortMatrix(GetNumRows(),GetNumColumns());
    for (int count=0; count<GetNumColumns(); count++)
        for (int count2=0; count2<GetNumRows(); count2++)
            shortMatrix(count2,count) = (short)mEigenMatrix(count2,count);

    return shortMatrix;
}

//! @brief ... converts any IntFullMatrix to an ShortFullMatrix
//! @return    converted matrix
template<>
FullMatrix<short> FullMatrix<int>::Convert2short()
{
    FullMatrix<short> shortMatrix(GetNumRows(),GetNumColumns());
    for (int count=0; count<GetNumColumns(); count++)
        for (int count2=0; count2<GetNumRows(); count2++)
            shortMatrix(count2,count) = (short)mEigenMatrix(count2,count);

    return shortMatrix;
}
//! @brief ... converts any ShortFullMatrix to a ShortFullMatrix
//! @return    converted matrix
template<>
FullMatrix<short> FullMatrix<short>::Convert2short()
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
void FullMatrix<short>::ImportFromSLangText(const char* fileName)
{
    throw MathException("[FullMatrix::importFromSLang] not implemented for this data-type.");
}

template<>
void FullMatrix<short>::ImportFromVtkASCIIFile(const char* rFileName)
{
    using namespace boost::spirit::classic;

    // open file
    std::ifstream file(rFileName, std::ios::in);
    if (file.is_open() == false)
    {
        throw MathException("[FullMatrix::ImportFromVtkASCIIFile] error opening file.");
    }

    // read header
    unsigned int imageDim[3];
    unsigned int numEntries(0);
    double imageSpacing[3], imageOrigin[3];
	assert(file.is_open());

	using namespace boost::spirit::classic;
	// read first four lines
	std::string line;
	getline (file, line);
	getline (file, line);
	getline (file, line);
	getline (file, line);
	// read dimension
	getline (file, line);
	if (parse(line.c_str(),("DIMENSIONS " >> uint_p[assign_a(imageDim[0])] >> ' '
								   >> uint_p[assign_a(imageDim[1])] >> ' '
								   >> uint_p[assign_a(imageDim[2])] >> *space_p)).full == false)
   {
		throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading dimension.");
   }
	// read spacing
	getline (file, line);
	if (parse(line.c_str(),("SPACING ">> real_p[assign_a(imageSpacing[0])] >> ' '
							>> real_p[assign_a(imageSpacing[1])] >> ' '
							>> real_p[assign_a(imageSpacing[2])] >> *space_p)).full == false)
	{
		throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading spacing.");
	}
	// read origin
	getline (file, line);
	if (parse(line.c_str(),("ORIGIN ">> real_p[assign_a(imageOrigin[0])] >> ' '
							 >> real_p[assign_a(imageOrigin[1])] >> ' '
							 >> real_p[assign_a(imageOrigin[2])] >> *space_p)).full == false)
	{
		 throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading origin.");
	}
	// read number of entries
	getline (file, line);
	if (parse(line.c_str(),("POINT_DATA ">> uint_p[assign_a(numEntries)] >>  *space_p)).full == false)
		   {
			   throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading number of entries.");
		   }
	std::cout<<"numEntries"<<numEntries<<"\n";
	// read data type
	getline (file, line);
	// read empty line
	getline (file, line);

    //this->ImportFromVtkASCIIFileReadHeader(file, imageDim, imageSpacing, imageOrigin, objectNumEntries);

    if (imageDim[2] = 0)
    {
        throw MathException("[FullMatrix::importFromVTkASCIIFile] only three dimension implemented.");
    }
    this->Resize(numEntries, 1);

    // read entries
    std::vector<unsigned short> imageValues;
  //  std::vector<unsigned short> shortImageValues;

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
    for (unsigned int count = 0; count<numEntries; count++)
    {
        this->mEigenMatrix(count, 0) = imageValues[count];
    }
    // close file
   file.close();
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
    unsigned int numEntries(0);

    std::string line;
    for (unsigned int count=0;count<7;count++)
    {
        getline (file, line);
    }

    // read number of entries
    getline (file, line);
    if (parse(line.c_str(),("POINT_DATA ">> uint_p[assign_a(numEntries)] >>  *space_p)).full == false)
           {
               throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading number of entries.");
           }
    // read data type
    getline (file, line);
    // read empty line
    getline (file, line);


    /*
    // read header
    unsigned int imageDim[3];
    double imageSpacing[3], imageOrigin[3];
	assert(file.is_open());

	// read first four lines
	std::string line;
	getline (file, line);
	getline (file, line);
	getline (file, line);
	getline (file, line);
	// read dimension
	getline (file, line);
	if (parse(line.c_str(),("DIMENSIONS " >> uint_p[assign_a(imageDim[0])] >> ' '
								   >> uint_p[assign_a(imageDim[1])] >> ' '
								   >> uint_p[assign_a(imageDim[2])] >> *space_p)).full == false)
   {
		throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading dimension.");
   }
	// read spacing
	getline (file, line);
	if (parse(line.c_str(),("SPACING ">> real_p[assign_a(imageSpacing[0])] >> ' '
							>> real_p[assign_a(imageSpacing[1])] >> ' '
							>> real_p[assign_a(imageSpacing[2])] >> *space_p)).full == false)
	{
		throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading spacing.");
	}
	// read origin
	getline (file, line);
	if (parse(line.c_str(),("ORIGIN ">> real_p[assign_a(imageOrigin[0])] >> ' '
							 >> real_p[assign_a(imageOrigin[1])] >> ' '
							 >> real_p[assign_a(imageOrigin[2])] >> *space_p)).full == false)
	{
		 throw MathException("[Matrix::importFromVtkASCIIFileReadHeader]error reading origin.");
	}
*/
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
    for (unsigned int count = 0; count<numEntries; count++)
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
FullMatrix<int> FullMatrix<int>::Abs()
{
	return FullMatrix<int> ( mEigenMatrix.cwise().abs() );
}

template<>
 FullMatrix<double> FullMatrix<double>::Abs()
{
	return FullMatrix<double> ( mEigenMatrix.cwise().abs() );
}

template<>
FullMatrix<short> FullMatrix<short>::Abs()
{
	throw MathException("[FullMatrix::Abs] not implemented in Eigen for short data-type.");
}
}
