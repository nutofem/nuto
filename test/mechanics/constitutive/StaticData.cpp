#define BOOST_TEST_MODULE StaticDataTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

//using namespace NuTo;
BOOST_AUTO_TEST_CASE(compenent)
{
}
BOOST_AUTO_TEST_CASE(composite)
{
}
BOOST_AUTO_TEST_CASE(primitive)
{
}
