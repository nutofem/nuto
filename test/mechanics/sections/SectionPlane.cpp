#include "BoostUnitTest.h"
#include <boost/test/output_test_stream.hpp>
#include <fstream>

#include "mechanics/MechanicsException.h"
#include "mechanics/sections/SectionPlane.h"

BOOST_AUTO_TEST_CASE(CreateAndPrintPlaneStrain)
{
    boost::test_tools::output_test_stream output;
    auto section = NuTo::SectionPlane::Create(42.0, true);
    output << *section;

    std::string expected = "    Plane section with thickness: 42\n    Section type is plane strain.\n";
    BOOST_CHECK(output.is_equal(expected));

    // plane section has no area; therefore it should throw
    BOOST_CHECK_THROW(section->GetArea(), NuTo::MechanicsException);
}

BOOST_AUTO_TEST_CASE(CreateAndPrintPlaneStress)
{
    boost::test_tools::output_test_stream output;
    auto section = NuTo::SectionPlane::Create(42.0, false);
    output << *section;

    std::string expected = "    Plane section with thickness: 42\n    Section type is plane stress.\n";
    BOOST_CHECK(output.is_equal(expected));
}
