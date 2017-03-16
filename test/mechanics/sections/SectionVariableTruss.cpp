#include "BoostUnitTest.h"
#include <boost/test/output_test_stream.hpp>
#include <fstream>

#include "mechanics/sections/SectionVariableTruss.h"

BOOST_AUTO_TEST_CASE(CreateAndPrint)
{
    boost::test_tools::output_test_stream output;
    auto section = NuTo::SectionVariableTruss::Create(0.0, 0.0, 0.0, 0.0);
    output << *section;

    std::string expected = "    Variable truss section.\n";
    BOOST_CHECK(output.is_equal(expected));
}

