//
// Created by Thomas Titscher on 10/21/16.
//
#define BOOST_TEST_MODULE DataContainerTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "nuto/mechanics/constitutive/staticData/DataContainer.h"
#include <type_traits>
// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}


namespace DataContainerTest
{
using namespace NuTo::Constitutive::StaticData;

//! @brief Test class that only needs to provide its StaticDataType
class TestLaw
{
public:
    typedef int StaticDataType;
};


BOOST_AUTO_TEST_CASE(DataContainer_Copy_Move)
{
    BOOST_CHECK(std::is_copy_constructible<DataContainer<TestLaw>>::value);
    BOOST_CHECK(std::is_move_constructible<DataContainer<TestLaw>>::value);

    BOOST_CHECK(std::is_copy_assignable<DataContainer<TestLaw>>::value);
    BOOST_CHECK(std::is_move_assignable<DataContainer<TestLaw>>::value);
}

BOOST_AUTO_TEST_CASE(DataContainerConstruction)
{
    {
        DataContainer<TestLaw> data(6174);
        BOOST_CHECK_EQUAL(data.GetNumData(), 1);
        BOOST_CHECK_EQUAL(data.GetData(), 6174);
    }
    {
        std::vector<typename TestLaw::StaticDataType> dataVector({0,1,2});
        DataContainer<TestLaw> data(dataVector);
        BOOST_CHECK_EQUAL(data.GetNumData(), 3);
        BOOST_CHECK_EQUAL(data.GetData(0), 0);
        BOOST_CHECK_EQUAL(data.GetData(1), 1);
        BOOST_CHECK_EQUAL(data.GetData(2), 2);
    }
    {
        DataContainer<TestLaw> data({0,1,2});
        BOOST_CHECK_EQUAL(data.GetNumData(), 3);
        BOOST_CHECK_EQUAL(data.GetData(0), 0);
        BOOST_CHECK_EQUAL(data.GetData(1), 1);
        BOOST_CHECK_EQUAL(data.GetData(2), 2);
    }
}

BOOST_AUTO_TEST_CASE(DataContainerShits)
{
    DataContainer<TestLaw> data(6174);

    // shifts not possible with only one set of static data
    BOOST_CHECK_THROW(data.ShiftToPast(), NuTo::MechanicsException);
    BOOST_CHECK_THROW(data.ShiftToFuture(), NuTo::MechanicsException);

    data.AllocateAdditionalData(4);
    BOOST_CHECK_EQUAL(data.GetData(4), 6174);
    BOOST_CHECK_EQUAL(data.GetNumData(), 5);

    data.SetData(1337);
    data.ShiftToPast(); // 1337 should be at [1]
    data.ShiftToPast(); // 1337 should be at [2]
    BOOST_CHECK_EQUAL(data.GetData(2), 1337);
    BOOST_CHECK_EQUAL(data.GetNumData(), 5); // should not alter number of static data

    data.ShiftToFuture(); // 1337 should be back at [1]
    data.ShiftToFuture(); // 1337 should be back at [0]
    BOOST_CHECK_EQUAL(data.GetData(), 1337);
    BOOST_CHECK_EQUAL(data.GetNumData(), 5);
}

BOOST_AUTO_TEST_CASE(DataContainerSerialze)
{
    DataContainer<TestLaw> data(3);
    data.AllocateAdditionalData(2);
    data.ShiftToFuture();
    data.SetData(2);
    data.ShiftToFuture();
    data.SetData(1);

    DataContainer<TestLaw> dataFromFile(0);
    dataFromFile.AllocateAdditionalData(2);

    {   // write
        NuTo::SerializeStreamOut s("DataContainerSerializeText.dat", false);
        s << data;
    }

    NuTo::SerializeStreamIn s("DataContainerSerializeText.dat", false);
    s >> dataFromFile;

    for (unsigned int i = 0; i < data.GetNumData(); ++i)
        BOOST_CHECK_EQUAL(dataFromFile.GetData(i),data.GetData(i));
}

} // namespace DataContainerTest