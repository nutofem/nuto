#include "BoostUnitTest.h"
#include <boost/test/output_test_stream.hpp>
#include <fstream>

#include "mechanics/MechanicsException.h"
#include "mechanics/sections/SectionTruss.h"

BOOST_AUTO_TEST_CASE(CreateAndPrint)
{
    boost::test_tools::output_test_stream output;
    auto section = NuTo::SectionTruss::Create(42.0);
    output << *section;

    std::string expected = "    Truss section with area 42\n";
    BOOST_CHECK(output.is_equal(expected));

    // truss has neither thickness or circumference, nor is it a plane, therefore it should throw
    BOOST_CHECK_THROW(section->GetThickness(), NuTo::MechanicsException);
    BOOST_CHECK_THROW(section->GetCircumference(), NuTo::MechanicsException);
    BOOST_CHECK_THROW(section->IsPlaneStrain(), NuTo::MechanicsException);
}

