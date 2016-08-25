// $Id$

#pragma once

#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <assert.h>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif// ENABLE_SERIALIZATION

#include <boost/spirit/include/classic_core.hpp>

#include "nuto/base/NuToObject.h"
#include "nuto/math/MathException.h"
#include "nuto/math/NuToMath.h"
#include "nuto/math/Operator.h"

namespace NuTo
{

//! @author Jörg F. Unger, ISM
//! @date July 2008
//! @brief ... standard abstract class for all matrices in NuTo
template <class T>
class Matrix : public NuToObject
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif// ENABLE_SERIALIZATION
public:
    Matrix() : NuToObject()
    {
    }

    //! @brief ... return the number of rows of the matrix
    //! @return string ... number of rows
    virtual int GetNumRows()const=0;

    //! @brief ... return the number of columns of the matrix
    //! @return string ... number of columns
    virtual int GetNumColumns()const=0;

    //! @brief ... add nonzero entry to matrix
    //! @param row ... row of the nonzero entry (zero based indexing!!!)
    //! @param column ... column of the nonzero entry (zero based indexing!!!)
    //! @param value ... value of the nonzero entry
    virtual void AddValue(int row, int column, const T& value) = 0;

    //! @brief ... resize the matrix and set all entries to zero
    //! @param row ... rows
    //! @param column ... columns
    virtual void Resize(int row, int column) = 0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject);
    }
#endif// ENABLE_SERIALIZATION

    //! @brief performs a monadic operator on all matrix entries
    //! @param rMOperator        Monadic Operator
    virtual void Map(const NuTo::MonadicOperator<T>* rMOperator)=0;

    //! @brief performs a dyadic operator on all matrix entries with another given value
    //! @param rDOperator        Dyadic Operator
    virtual void Map(const NuTo::DyadicOperator<T>* rDOperator, const T& rValue)=0;

    //! @brief ... convert the entries of the matrix to strings (for file export and std::cout)
    //! @param value ... entry of the matrix
    //! @return string ... representation of the entry
    static std::string Convert2String(T value);

    //! @brief ... convert the entries of the matrix to strings (for file export and std::cout)
    //! @param value ... entry of the matrix
    //! @param scientific ... scientific or fixed point notation
    //! @param precision ... number of digits after the comma
    //! @param width ... total width
    //! @return string ... representation of the entry
    std::string Convert2String(T value, bool scientific, int precision, int width) const;

    //! @brief ... convert a string to an entries of the matrix (for file import)
    //! @param s ...  string to convert
    //! @return ... extracted value
    static T ConvertFromString(const std::string s)
    {
        T toT;
        std::istringstream(s)>>toT;
        return toT;
    }

    //! @brief ... convert an integer to strings (for file export and std::cout)
    //! @param i ... integer value
    //! @return string ... representation of the integer
    static std::string Int2String(int i)
    {
        std::ostringstream format_message;
        format_message << i;
        return format_message.str();
    }

    //! @brief ... convert a double to strings (for file export and std::cout)
    //! @param value ... double value
    //! @return string ... representation of the double
    static std::string Double2String(double value)
    {
        std::ostringstream format_message;
        format_message << std::setprecision(12) << std::setw(15) << std::showpoint << value;
        return format_message.str();
    }

protected:

    //! @brief ... resizes and fills the matrix rMatrix with rNumValues random values
    //! @param rMatrix ... Matrix<T>
    //! @param rDensity ... approximate density = numValues / (rNumRows*rNumColumns)
    //! @param rSeed ... random seed
    static void FillMatrixRandom(Matrix& rMatrix, double rDensity, int rSeed)
    {
        std::mt19937 gen(rSeed);
        std::uniform_real_distribution<double> value_distribution(0., 10.);
        std::uniform_int_distribution<int> row_distribution(0, rMatrix.GetNumRows() - 1);
        std::uniform_int_distribution<int> col_distribution(0, rMatrix.GetNumColumns() - 1);

        int numValues = (int) (rDensity * rMatrix.GetNumRows() * rMatrix.GetNumColumns());

        for (int i = 0; i < numValues; ++i)
        {
            int row = row_distribution(gen);
            int col = col_distribution(gen);
            double val = value_distribution(gen);
            rMatrix.AddValue(row, col, val);
        }
    }

    enum SLangObjectKind
    {
        SLANG_NOKIND=0,
        SLANG_SCALAR,
        SLANG_VECTOR,
        SLANG_MATRIX,
        SLANG_UPPER_TRIANGLE,
        SLANG_LOWER_TRIANGLE,
        SLANG_SKYLINE_MATRIX,
        SLANG_COMPACT_MATRIX,
        SLANG_GENERAL_COMPACT_MATRIX,
        SLANG_RAGGED_MATRIX
    };
    enum SLangObjectType
    {
        SLANG_NOTYPE=0,
        SLANG_INTEGER,
        SLANG_REAL,
        SLANG_CHARACTER,
        SLANG_COMPLEX,
        SLANG_SINGLE_TYPE
    };

    void ImportFromSLangTextReadHeader(std::ifstream& file, unsigned int& SLangVersion, unsigned int& type, unsigned int& kind, unsigned int& first, unsigned int& second, unsigned int& numEntries)
    {
        assert(file.is_open());

        using namespace boost::spirit::classic;
        // read first line SLang version
        std::string line;
        getline (file, line);
        unsigned int SLangVersionMajor(0), SLangVersionMinor(0), SLangVersionDevelop(0);
        if (parse(line.c_str(),("SLang " >> uint_p[assign_a(SLangVersionMajor)] >> '.'
                                >> uint_p[assign_a(SLangVersionMinor)] >> '.'
                                >> uint_p[assign_a(SLangVersionDevelop)] >> *space_p)).full == false)
        {
            throw MathException("[Matrix::importFromSLangTextReadHeader]error reading SLang version.");
        }
        SLangVersion = SLangVersionMajor * 100 + SLangVersionMinor * 10 + SLangVersionDevelop;

        // read next two lines
        getline (file, line);
        getline (file, line);

        // read object info
        getline (file, line);
        if (parse(line.c_str(),("Object info: " >> uint_p[assign_a(type)] >> ' '
                                >> uint_p[assign_a(kind)] >> ' '
                                >> uint_p[assign_a(first)] >> ' '
                                >> uint_p[assign_a(second)] >> ' '
                                >> uint_p[assign_a(numEntries)] >> *space_p)).full == false)
        {
            throw MathException("[Matrix::importFromSLangTextReadHeader]error reading object info.");
        }

        // read empty line
        getline (file, line);
    }
    /*
    	virtual void Max(T& rResultOutput)=0;
    	virtual void Max(int& rRowOutput, int& rColumnOutput, T& rResultOutput)=0;
    	virtual void Min(T& rResultOutput)=0;
    	virtual void Min(int& rRowOutput, int& rColumnOutput, T& rResultOutput)=0;
    */
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Matrix<int>)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Matrix<double>)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
