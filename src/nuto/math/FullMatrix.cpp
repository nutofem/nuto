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

}
