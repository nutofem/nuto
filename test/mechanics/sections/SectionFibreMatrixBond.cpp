#include "BoostUnitTest.h"
#include <boost/test/output_test_stream.hpp>
#include <fstream>

#include "mechanics/sections/SectionFibreMatrixBond.h"

BOOST_AUTO_TEST_CASE(CreateAndPrint)
{
    boost::test_tools::output_test_stream output;
    auto section = NuTo::SectionFibreMatrixBond::Create(42.0);
    output << *section;

    std::string expected = "    Fibre matrix bond section with circumference: 42\n";
    BOOST_CHECK(output.is_equal(expected));
}

