#define BOOST_TEST_MODULE SerializeDataTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <memory>
#include <eigen3/Eigen/Dense>

#include "nuto/base/Timer.h"

#include "nuto/math/FullVector.h"

#include "nuto/base/serializeStream/SerializeStreamOut.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"


// needed for building with clang when boost test has been built with gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size()));
}

namespace SerializeDataTest
{

template <typename T>
void WriteAndRead(const std::string& rFile, bool rIsBinary, const T& rValue, T& rValueFromFile)
{
    // write
    {
        NuTo::SerializeStreamOut streamOut(rFile, rIsBinary);
        streamOut.WriteData(rValue);
    } // out stream out of scope and f***ing closes the f***ing file. Damnit.

    // read
    NuTo::SerializeStreamIn streamIn(rFile, rIsBinary);
    streamIn.ReadData(rValueFromFile);
}

void Double(const std::string& rFile, bool rIsBinary)
{
    double value = 1337.6174e42;
    double valueFromFile = -42;

    WriteAndRead(rFile, rIsBinary, value, valueFromFile);
    BOOST_CHECK_CLOSE(value, valueFromFile, 1.e-10);
}

void EigenMatrixDynamic(const std::string& rFile, bool rIsBinary)
{
    Eigen::MatrixXd values = Eigen::MatrixXd::Random(3, 15);
    Eigen::MatrixXd valuesFromFile(3, 15);

    WriteAndRead(rFile, rIsBinary, values, valuesFromFile);
    BOOST_CHECK_CLOSE(values.norm(), valuesFromFile.norm(), 1.e-10);
}

void EigenMatrixFixed(const std::string& rFile, bool rIsBinary)
{
    Eigen::Matrix3d values = Eigen::Matrix3d::Random();
    Eigen::Matrix3d valuesFromFile;

    WriteAndRead(rFile, rIsBinary, values, valuesFromFile);
    BOOST_CHECK_CLOSE(values.norm(), valuesFromFile.norm(), 1.e-10);
}

void FullVector(const std::string& rFile, bool rIsBinary)
{
    NuTo::FullVector<double> values = NuTo::FullVector<double>::Random(4);
    NuTo::FullVector<double> valuesFromFile(4);

    WriteAndRead(rFile, rIsBinary, values, valuesFromFile);
    BOOST_CHECK_CLOSE(values.norm(), valuesFromFile.norm(), 1.e-10);
}

class CompoundDataBase
{
public:
    virtual ~CompoundDataBase() = default;
    friend NuTo::SerializeStreamOut& operator<<(NuTo::SerializeStreamOut& rStream, CompoundDataBase& rData)
    {
        rData.WriteMyData(rStream);
        return rStream;
    }

    friend NuTo::SerializeStreamIn& operator>>(NuTo::SerializeStreamIn& rStream, CompoundDataBase& rData)
    {
        rData.ReadMyData(rStream);
        return rStream;
    }

    virtual void WriteMyData(NuTo::SerializeStreamOut& rStream) = 0;
    virtual void ReadMyData(NuTo::SerializeStreamIn& rStream) = 0;
};

class CompoundData : public CompoundDataBase
{
public:

    double mScalar;
    Eigen::Vector3d mVector;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mMatrix;

    void SetRandom(int rRows, int rCols)
    {
        mScalar = 42.1337; // Hey... this number is PRETTY random... Have you thought of it?
        mVector = Eigen::Vector3d::Random();
        mMatrix = Eigen::MatrixXd::Random(rRows, rCols);
    }

    void SetZero(int rRows = 0, int rCols = 0)
    {
        mScalar = 0.;
        mVector.setZero();
        mMatrix.setZero(rRows, rCols);
    }

    void ReadMyData(NuTo::SerializeStreamIn& rStream) override
    {
        SerializeMyData(rStream);
    }
    void WriteMyData(NuTo::SerializeStreamOut& rStream) override
    {
        SerializeMyData(rStream);
    }

private:
    template <typename TStream>
    void SerializeMyData(TStream& rStream)
    {
        rStream.SerializeData(mScalar);
        rStream.SerializeData(mVector);
        rStream.SerializeData(mMatrix);
    }
};


void CheckCompoundData(const std::string& rFile, bool rIsBinary)
{
    CompoundData values1;
    CompoundData values2;

    values1.SetRandom(5,2);
    values2.SetRandom(3,7);

    CompoundData valuesFromFile1 = values1;
    CompoundData valuesFromFile2 = values2;

    valuesFromFile1.SetZero(5,2);
    valuesFromFile2.SetZero(3,7);

    // write
    {
        NuTo::SerializeStreamOut streamOut(rFile, rIsBinary);
        streamOut << values1;
        streamOut << values2;
    }

    // read
    NuTo::SerializeStreamIn streamIn(rFile, rIsBinary);
    streamIn >> valuesFromFile1;
    streamIn >> valuesFromFile2;

    BOOST_CHECK_CLOSE(values1.mVector.norm(),   valuesFromFile1.mVector.norm(), 1.e-10);
    BOOST_CHECK_CLOSE(values1.mMatrix.norm(),   valuesFromFile1.mMatrix.norm(), 1.e-10);

    BOOST_CHECK_CLOSE(values2.mVector.norm(),   valuesFromFile2.mVector.norm(), 1.e-10);
    BOOST_CHECK_CLOSE(values2.mMatrix.norm(),   valuesFromFile2.mMatrix.norm(), 1.e-10);
}

void CheckVectorCompoundData(const std::string& rFile, bool rIsBinary)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__);
    size_t num = 3;
    std::vector<std::unique_ptr<CompoundDataBase>> values(num);
    std::vector<std::unique_ptr<CompoundDataBase>> valuesFromFile(num);
    for (auto& value : values)
    {
        value = std::make_unique<CompoundData>();
        static_cast<CompoundData*>(&(*value))->SetRandom(5,3);
    }

    for (auto& valueFromFile : valuesFromFile)
    {
        valueFromFile = std::make_unique<CompoundData>();
        static_cast<CompoundData*>(&(*valueFromFile))->SetZero(5,3);
    }

    timer.Reset("Write");
    // write
    {
        NuTo::SerializeStreamOut streamOut(rFile, rIsBinary);
        for (auto& value : values)
            streamOut << *value;
    }
    timer.Reset("Read");
    // read
    NuTo::SerializeStreamIn streamIn(rFile, rIsBinary);
    for (auto& valueFromFile : valuesFromFile)
        streamIn >> *valueFromFile;

    timer.Reset("Cleanup");
    for (size_t i = 0; i < num; ++i)
    {
        BOOST_CHECK_CLOSE(static_cast<CompoundData*>(&*values[i])->mVector.norm(),
                          static_cast<CompoundData*>(&*valuesFromFile[i])->mVector.norm(), 1.e-10);
        BOOST_CHECK_CLOSE(static_cast<CompoundData*>(&*values[i])->mMatrix.norm(),
                          static_cast<CompoundData*>(&*valuesFromFile[i])->mMatrix.norm(), 1.e-10);
    }

}


}

BOOST_AUTO_TEST_CASE(SerializeDoubleText)               { SerializeDataTest::Double                 ("DoubleText.dat"                   , false);}
BOOST_AUTO_TEST_CASE(SerializeDoubleBinary)             { SerializeDataTest::Double                 ("DoubleBinary.dat"                 , true );}
BOOST_AUTO_TEST_CASE(SerializeEigenMatrixDynamicText)   { SerializeDataTest::EigenMatrixDynamic     ("EigenMatrixDynamicText.dat"       , false);}
BOOST_AUTO_TEST_CASE(SerializeEigenMatrixDynamicBinary) { SerializeDataTest::EigenMatrixDynamic     ("EigenMatrixDynamicBinary.dat"     , true );}
BOOST_AUTO_TEST_CASE(SerializeEigenMatrixFixedText)     { SerializeDataTest::EigenMatrixFixed       ("EigenMatrixFixedText.dat"         , false);}
BOOST_AUTO_TEST_CASE(SerializeEigenMatrixFixedBinary)   { SerializeDataTest::EigenMatrixFixed       ("EigenMatrixFixedBinary.dat"       , true );}
BOOST_AUTO_TEST_CASE(SerializeFullVectorText)           { SerializeDataTest::FullVector             ("FullVectorText.dat"               , false);}
BOOST_AUTO_TEST_CASE(SerializeFullVectorBinary)         { SerializeDataTest::FullVector             ("FullVectorBinary.dat"             , true );}
BOOST_AUTO_TEST_CASE(SerializeCompoundDataText)         { SerializeDataTest::CheckCompoundData      ("CompoundDataText.dat"             , false);}
BOOST_AUTO_TEST_CASE(SerializeCompoundDataBinary)       { SerializeDataTest::CheckCompoundData      ("CompoundDataBinary.dat"           , true );}
BOOST_AUTO_TEST_CASE(SerializeVectorCompoundDataText)   { SerializeDataTest::CheckVectorCompoundData("VectorCompoundDataText.dat"       , false);}
BOOST_AUTO_TEST_CASE(SerializeVectorCompoundDataBinary) { SerializeDataTest::CheckVectorCompoundData("VectorCompoundDataBinary.dat"     , true );}

