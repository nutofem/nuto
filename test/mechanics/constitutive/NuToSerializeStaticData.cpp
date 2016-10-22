#define BOOST_TEST_MODULE NuToSerializeStaticDataTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <eigen3/Eigen/Core>

#include "nuto/base/serializeStream/SerializeStreamOut.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"

// needed for building with clang when boost test has been built with gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size()));
}


namespace NuToSerializeStaticData
{

class A
{
public:

    A(const Eigen::Matrix2d& rData) : mData(rData) {}

    void NuToSerializeSave(NuTo::SerializeStreamOut& rStream)
    {
        rStream.Serialize(mData);
    }

    void NuToSerializeLoad(NuTo::SerializeStreamIn& rStream)
    {
        rStream.Serialize(mData);
    }

    Eigen::Matrix2d mData;

};

BOOST_AUTO_TEST_CASE(SerializeCompositeBase)
{
#warning "ToDo!"
}

} // namespace




