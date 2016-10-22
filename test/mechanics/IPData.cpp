//
// Created by Thomas Titscher on 10/22/16.
//

#define BOOST_TEST_MODULE IPDataTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLaw.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"


// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

using namespace NuTo::Constitutive;

class MockLaw
{
public:
    template<int TDim>
    NuTo::eError Evaluate(NuTo::ConstitutiveInputMap, NuTo::ConstitutiveOutputMap)
    {
        return NuTo::eError::SUCCESSFUL;
    }
};

BOOST_AUTO_TEST_CASE(withoutDataTest)
{
//    MockLaw law;
//    IPConstitutiveLawWithoutData<MockLaw> ip(law);
//    BOOST_CHECK(law == ip.GetConstitutiveLaw());

//    NuTo::ConstitutiveInputMap input;
//    NuTo::ConstitutiveOutputMap output;
//    ip.template Evaluate<1>(input, output);
}

