#define BOOST_TEST_MODULE IPConstitutiveLawTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawWithoutData.h"

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

using namespace NuTo::Constitutive;


BOOST_AUTO_TEST_CASE(IPConstitutiveLaw_Test)
{
//    TestLaw law;
//    auto iplaw = law.CreateIPLaw();
//
//
//    BOOST_CHECK(law == dynamic_cast<TestLaw&>(iplaw->GetConstitutiveLaw()));
//
//    NuTo::ConstitutiveInputMap input;
//    NuTo::ConstitutiveOutputMap output;
//    iplaw->Evaluate<1>(input, output);
}
