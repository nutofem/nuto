#include "BoostUnitTest.h"
#include "TypeTraits.h"

#include "mechanics/constitutive/staticData/DataContainer.h"

namespace DataContainerTest
{
using namespace NuTo::Constitutive::StaticData;

BOOST_AUTO_TEST_CASE(DataContainer_Copy_Move)
{
    NuTo::Test::Copy<DataContainer<int>>();
    NuTo::Test::Move<DataContainer<int>>();
}

BOOST_AUTO_TEST_CASE(DataContainerConstruction)
{
    {
        DataContainer<int> data(6174);
        BOOST_CHECK_EQUAL(data.GetNumData(), 1);
        BOOST_CHECK_EQUAL(data.GetData(), 6174);
    }
    {
        std::vector<int> dataVector({0, 1, 2});
        DataContainer<int> data(dataVector);
        BOOST_CHECK_EQUAL(data.GetNumData(), 3);
        BOOST_CHECK_EQUAL(data.GetData(0), 0);
        BOOST_CHECK_EQUAL(data.GetData(1), 1);
        BOOST_CHECK_EQUAL(data.GetData(2), 2);
    }
    {
        DataContainer<int> data({0, 1, 2});
        BOOST_CHECK_EQUAL(data.GetNumData(), 3);
        BOOST_CHECK_EQUAL(data.GetData(0), 0);
        BOOST_CHECK_EQUAL(data.GetData(1), 1);
        BOOST_CHECK_EQUAL(data.GetData(2), 2);
    }
}

BOOST_AUTO_TEST_CASE(DataContainerShits)
{
    DataContainer<int> data(6174);

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
    DataContainer<int> data(3);
    data.AllocateAdditionalData(2);
    data.ShiftToFuture();
    data.SetData(2);
    data.ShiftToFuture();
    data.SetData(1);

    DataContainer<int> dataFromFile(0);
    dataFromFile.AllocateAdditionalData(2);

    { // write
        NuTo::SerializeStreamOut s("DataContainerSerializeText.dat", false);
        s << data;
    }

    NuTo::SerializeStreamIn s("DataContainerSerializeText.dat", false);
    s >> dataFromFile;

    for (unsigned int i = 0; i < data.GetNumData(); ++i)
        BOOST_CHECK_EQUAL(dataFromFile.GetData(i), data.GetData(i));
}

} // namespace DataContainerTest
