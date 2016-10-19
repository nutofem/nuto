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

//! @remark provide rZero for proper initialization of the primitive types
template <typename T>
T WriteReadNumber(const std::string &rFile, bool rIsBinary, T &rValue, T rZero)
{
    // write
    {
        NuTo::SerializeStreamOut streamOut(rFile, rIsBinary);
        streamOut.Serialize(rValue);
    } // out stream out of scope and f***ing closes the f***ing file. Damnit.

    // read
    NuTo::SerializeStreamIn streamIn(rFile, rIsBinary);
    streamIn.Serialize(rZero);

    return rZero;
}

BOOST_AUTO_TEST_CASE(SerializeDouble)
{
    double d = 13.37e42;
    BOOST_CHECK_CLOSE(d, WriteReadNumber("DoubleText.dat", false, d, 0.), 1.e-10);
    BOOST_CHECK_CLOSE(d, WriteReadNumber("DoubleBinary.dat", true, d, 0.), 1.e-10);
}

template <typename T>
T WriteReadMatrix(const std::string &rFile, bool rIsBinary, T &rValue)
{
    // write
    {
        NuTo::SerializeStreamOut streamOut(rFile, rIsBinary);
        streamOut.Serialize(rValue);
    } // out stream out of scope and f***ing closes the f***ing file. Damnit.

    // read
    T valueFromFile = rValue;
    valueFromFile.setZero();
    NuTo::SerializeStreamIn streamIn(rFile, rIsBinary);
    streamIn.Serialize(valueFromFile);

    return valueFromFile;
}


BOOST_AUTO_TEST_CASE(SerializeEigen)
{
    {
        Eigen::MatrixXd m = Eigen::MatrixXd::Random(3, 2);
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("EigenDynamicText.dat", false, m).norm(), 1.e-10);
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("EigenDynamicBinary.dat", true, m).norm(), 1.e-10);
    }
    {
        Eigen::Vector3d m = Eigen::Vector3d::Random();
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("EigenFixedText.dat", false, m).norm(), 1.e-10);
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("EigenFixedBinary.dat", true, m).norm(), 1.e-10);
    }
}
BOOST_AUTO_TEST_CASE(SerializeNuTo)
{
    {
        NuTo::FullVector<double, Eigen::Dynamic> m = Eigen::VectorXd::Random(2);
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("NuToVectorDynamicText.dat", false, m).norm(), 1.e-10);
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("NuToVectorDynamicBinary.dat", true, m).norm(), 1.e-10);
    }
    {
        NuTo::FullMatrix<double, 2,2> m = Eigen::Matrix2d::Random();
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("NuToMatrixFixedText.dat", false, m).norm(), 1.e-10);
        BOOST_CHECK_CLOSE(m.norm(), WriteReadMatrix("NuToMatrixFixedBinary.dat", true, m).norm(), 1.e-10);
    }
}



class CompoundDataBase
{
public:
    virtual ~CompoundDataBase() = default;

    virtual void NuToSerializeSave(NuTo::SerializeStreamOut& rStream)
    {
        rStream.Serialize(mScalar);
    }
    virtual void NuToSerializeLoad(NuTo::SerializeStreamIn& rStream)
    {
        rStream.Serialize(mScalar);
    }

    double mScalar;
};

class CompoundData : public CompoundDataBase
{
public:

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

    virtual void NuToSerializeSave(NuTo::SerializeStreamOut& rStream) override
    {
        CompoundDataBase::NuToSerializeSave(rStream);  // explicitly call the base class serialization
        SerializeMe(rStream);
    }
    virtual void NuToSerializeLoad(NuTo::SerializeStreamIn& rStream) override
    {
        CompoundDataBase::NuToSerializeLoad(rStream);  // explicitly call the base class serialization
        SerializeMe(rStream);
    }

private:
    template <typename TStream>
    void SerializeMe(TStream& rStream)
    {
        rStream.Serialize(mVector);
        rStream.Serialize(mMatrix);
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

    BOOST_CHECK_CLOSE(values1.mScalar,          valuesFromFile1.mScalar,        1.e-10);
    BOOST_CHECK_CLOSE(values1.mVector.norm(),   valuesFromFile1.mVector.norm(), 1.e-10);
    BOOST_CHECK_CLOSE(values1.mMatrix.norm(),   valuesFromFile1.mMatrix.norm(), 1.e-10);

    BOOST_CHECK_CLOSE(values2.mScalar,          valuesFromFile2.mScalar,        1.e-10);
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
        BOOST_CHECK_CLOSE(static_cast<CompoundData*>(&*values[i])->mScalar,
                          static_cast<CompoundData*>(&*valuesFromFile[i])->mScalar, 1.e-10);
        BOOST_CHECK_CLOSE(static_cast<CompoundData*>(&*values[i])->mVector.norm(),
                          static_cast<CompoundData*>(&*valuesFromFile[i])->mVector.norm(), 1.e-10);
        BOOST_CHECK_CLOSE(static_cast<CompoundData*>(&*values[i])->mMatrix.norm(),
                          static_cast<CompoundData*>(&*valuesFromFile[i])->mMatrix.norm(), 1.e-10);
    }

}

BOOST_AUTO_TEST_CASE(SerializeCompoundDataText)         { CheckCompoundData      ("CompoundDataText.dat"             , false);}
BOOST_AUTO_TEST_CASE(SerializeCompoundDataBinary)       { CheckCompoundData      ("CompoundDataBinary.dat"           , true );}
BOOST_AUTO_TEST_CASE(SerializeVectorCompoundDataText)   { CheckVectorCompoundData("VectorCompoundDataText.dat"       , false);}
BOOST_AUTO_TEST_CASE(SerializeVectorCompoundDataBinary) { CheckVectorCompoundData("VectorCompoundDataBinary.dat"     , true );}

}
