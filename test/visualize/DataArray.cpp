#include "BoostUnitTest.h"
#include "nuto/visualize/DataArray.h"
#include <sstream>
#include <fstream>
#include <iostream>

BOOST_AUTO_TEST_CASE(DataArrayFormatString)
{
    NuTo::Visualize::DataArray<double> aDouble("", 3, {});
    NuTo::Visualize::DataArray<float> aFloat("", 3, {});
    NuTo::Visualize::DataArray<unsigned> aUnsigned("", 3, {});
    NuTo::Visualize::DataArray<uint8_t> aShort("", 3, {});
    BOOST_CHECK_EQUAL(aDouble.VtkTypeString(), "Float64");
    BOOST_CHECK_EQUAL(aFloat.VtkTypeString(), "Float32");
    BOOST_CHECK_EQUAL(aUnsigned.VtkTypeString(), "UInt32");
    BOOST_CHECK_EQUAL(aShort.VtkTypeString(), "UInt8");
}

BOOST_AUTO_TEST_CASE(DataArrayOffset)
{
    NuTo::Visualize::DataArray<double> aDouble("", 3, {0, 0, 0, 0, 0, 0});
    NuTo::Visualize::DataArray<uint8_t> aShort("", 2, {0, 0, 0, 0});

    // header + #data * sizeof(double)
    BOOST_CHECK_EQUAL(aDouble.Offset(), 8 + 6 * 8);
    BOOST_CHECK_EQUAL(aShort.Offset(), 8 + 4 * 1);
}

BOOST_AUTO_TEST_CASE(DataArrayCommonHeader)
{
    NuTo::Visualize::DataArray<unsigned> dataArray("Name", 3, {1, 2, 3, 4, 5, 6});
    std::stringstream ss;
    dataArray.WriteCommonHeader(ss);
    std::string header = ss.str();

    BOOST_CHECK(header.find(R"(<DataArray)") != std::string::npos);
    BOOST_CHECK(header.find(R"( type="UInt32" )") != std::string::npos);
    BOOST_CHECK(header.find(R"( Name="Name" )") != std::string::npos);
    BOOST_CHECK(header.find(R"( NumberOfComponents="3" )") != std::string::npos);

    NuTo::Visualize::DataArray<unsigned> dataArray2("Name", 0, {1, 2, 3, 4, 5, 6});
    ss.clear();
    dataArray2.WriteCommonHeader(ss);
    // zero components: should not be in the header
    BOOST_CHECK(ss.str().find(R"( NumberOfComponents="0" )") == std::string::npos);
}


BOOST_AUTO_TEST_CASE(DataArrayAscii)
{
    NuTo::Visualize::DataArray<unsigned> dataArray("Name", 3, {1, 2, 3, 4, 5, 6});
    std::stringstream ss;
    dataArray.WriteAscii(ss);
    std::string asciiData = ss.str();

    BOOST_CHECK(asciiData.find(R"( format="ascii" )") != std::string::npos);
    BOOST_CHECK(asciiData.find(R"( 1 2 3 4 5 6)") != std::string::npos);
    BOOST_CHECK(asciiData.find(R"(</DataArray>)") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(DataArrayAsciiSmallDouble)
{
    NuTo::Visualize::DataArray<double> dataArray("Name", 1, {5.1e-8});
    std::stringstream ss;
    dataArray.WriteAscii(ss);
    std::string asciiData = ss.str();

    BOOST_CHECK(asciiData.find(R"(5.1e-08)") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(DataArrayAsciiUint8t)
{
    uint8_t val = 3;
    NuTo::Visualize::DataArray<uint8_t> dataArray("Name", 1, {val});
    std::stringstream ss;
    dataArray.WriteAscii(ss);
    std::string asciiData = ss.str();

    BOOST_CHECK(asciiData.find(R"(3)") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(DataArrayBinaryHeader)
{
    NuTo::Visualize::DataArray<unsigned> dataArray("Name", 3, {1, 2, 3, 4, 5, 6});
    std::stringstream ss;
    int offset = 12;
    dataArray.WriteBinaryHeader(ss, &offset);
    BOOST_CHECK_EQUAL(offset, 12 + 8 + 6 * sizeof(unsigned));

    std::string header = ss.str();

    BOOST_CHECK(header.find(R"( format="appended" )") != std::string::npos);
    BOOST_CHECK(header.find(R"( offset="12" )") != std::string::npos);
    BOOST_CHECK(header.find(R"(/>)") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(DataArrayBinary)
{
    NuTo::Visualize::DataArray<double> dataArray("Name", 3, {1, 2, 3, 4, 5, 6});
    std::ofstream file("DataArrayBinary.dat");
    dataArray.WriteBinaryData(file);
    file.close();

    std::ifstream binaryFile("DataArrayBinary.dat");
    uint64_t size;
    binaryFile.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    BOOST_CHECK_EQUAL(size, 6 * sizeof(double));

    std::vector<double> data(6);
    binaryFile.read(reinterpret_cast<char*>(data.data()), size);
    BoostUnitTest::CheckVector(data, std::vector<double>({1, 2, 3, 4, 5, 6}), 6);
}
