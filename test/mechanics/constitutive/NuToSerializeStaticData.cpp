#define BOOST_TEST_MODULE NuToSerializeStaticDataTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "nuto/base/Timer.h"

#include "nuto/base/serializeStream/SerializeStreamOut.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"

#include "nuto/mechanics/constitutive/staticData/Leaf.h"

// needed for building with clang when boost test has been built with gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size()));
}


namespace NuToSerializeStaticData
{

template <typename T>
void SerializeLeaf(T value)
{
//    auto leaf = NuTo::Constitutive::StaticData::Leaf<T>::Create(value);
//    std::cout << leaf->GetData();
}


} // namespace




BOOST_AUTO_TEST_CASE(NuToSerializeLeafDouble) {NuToSerializeStaticData::SerializeLeaf<double>(5.);}