//
// Created by ttitsche on 1/26/17.
//

#define BOOST_TEST_MODULE MeshCompanionTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"


BOOST_AUTO_TEST_CASE(ImportBinaryFromGmsh)
{
    NuTo::Structure s(3);
    NuTo::MeshCompanion::ImportFromGmsh(s, "MeshCompanionGmsh.msh");
}