#include "nuto/mechanics/constitutive/staticData/Leaf.h"
#include "nuto/mechanics/constitutive/staticData/DataMoistureTransport.h"

#define BOOST_TEST_MODULE LeafTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

#include <iostream>

using namespace NuTo::Constitutive::StaticData;

BOOST_AUTO_TEST_CASE(add_and_retrieve_scalars)
{
    Leaf<double> historyVariable;

    historyVariable.SetData(42); // 42 is now at the "current" static data
    historyVariable.AllocateAdditionalData(4);
    BOOST_CHECK_EQUAL(historyVariable.GetData(4), 42);
    BOOST_CHECK_EQUAL(historyVariable.GetNumData(), 5);

    historyVariable.SetData(1337);
    historyVariable.ShiftToPast(); // 1337 should be at [1]
    historyVariable.ShiftToPast(); // 1337 should be at [2]
    BOOST_CHECK_EQUAL(historyVariable.GetData(2), 1337);
    BOOST_CHECK_EQUAL(historyVariable.GetNumData(), 5); // should not alter number of static data

    historyVariable.ShiftToFuture(); // 1337 should be back at [1]
    historyVariable.ShiftToFuture(); // 1337 should be back at [0]
    BOOST_CHECK_EQUAL(historyVariable.GetData(), 1337);
    BOOST_CHECK_EQUAL(historyVariable.GetNumData(), 5);

    Leaf<double>* historyVariableTwo = historyVariable.Clone();
    BOOST_CHECK(historyVariable == *historyVariableTwo);

    historyVariableTwo->SetData(27);
    BOOST_CHECK(historyVariable != *historyVariableTwo);
}


BOOST_AUTO_TEST_CASE(add_and_retrieve_custom_types)
{
    Leaf<DataMoistureTransport> moistureData;

    moistureData.SetData(DataMoistureTransport());
    BOOST_CHECK_EQUAL(moistureData.GetNumData(), 1);

    auto currentData = moistureData.GetData();
    currentData.SetLastRelHumValue(7.0);
    moistureData.SetData(currentData);
    BOOST_CHECK_EQUAL(moistureData.GetData(0).GetLastRelHumValue(), 7.0);

    moistureData.AllocateAdditionalData(3);
    BOOST_CHECK_EQUAL(moistureData.GetNumData(), 4);

    moistureData.GetData(3).SetDesorption(false);
    BOOST_CHECK(!moistureData.GetData(3).IsDesorption());

    moistureData.ShiftToPast();
    moistureData.ShiftToPast();
    moistureData.ShiftToFuture();
    moistureData.ShiftToFuture();

    Leaf<DataMoistureTransport>* moistureDataTwo = moistureData.Clone();

    BOOST_CHECK(moistureData == *moistureDataTwo);

    currentData.SetLastRelHumValue(500.4);
    moistureDataTwo->SetData(currentData);
    BOOST_CHECK(moistureData != *moistureDataTwo);
}
