#pragma once
/*
 * collect boost::unittest boilerplate here
 */

#define BOOST_TEST_MODULE DummyTestName
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}